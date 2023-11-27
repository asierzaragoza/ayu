# This requires the prediction_results.txt file
import pandas as pd
import numpy as np
import os
import itertools
import pickle
import scipy.special
import skbio.stats.composition as skbio_comp
import ayu.preprocessing
from Bio.SeqIO.FastaIO import SimpleFastaParser

#TMBed
def parse_tmbed_file(tmbed_file):
    '''PRD: Predicted type, the other keys are the prediction values for each type'''
    tmbed_master_list = []
    with open(tmbed_file) as in_handle:
        selected_dict = {'prot_ID':None}
        for line in in_handle:
            if line[0] == '>':
                if selected_dict['prot_ID'] is not None:
                    tmbed_master_list.append(selected_dict)
                selected_dict = {'prot_ID':line.rstrip('\n')[1:],
                    'PRD':[], 'P_B':[], 'P_H':[], 'P_S':[], 'P_i':[],
                    'P_o':[]}
                in_handle.readline()

            else:
                splitLine = line.rstrip('\n').split('\t')
                selected_dict['PRD'].append(splitLine[1])
                selected_dict['P_B'].append(float(splitLine[2]))
                selected_dict['P_H'].append(float(splitLine[3]))
                selected_dict['P_S'].append(float(splitLine[4]))
                selected_dict['P_i'].append(float(splitLine[5]))
                selected_dict['P_o'].append(float(splitLine[6]))
    
    return tmbed_master_list

def get_avg_values(tmbed_list):
    result_list = []
    for tmbed_dict in tmbed_list:
        tmbed_avg_dict = {'prot_ID':tmbed_dict['prot_ID']}
        for key in tmbed_dict:
            if key not in set(['prot_ID', 'PRD']):
                tmbed_avg_dict[key] = round(sum(tmbed_dict[key]) / len(tmbed_dict[key]), 4)
        result_list.append(tmbed_avg_dict)
    return result_list

def add_slr_to_tmbed(tmbed_file, out_file):
    def closure_func(x):
        return pd.Series(skbio_comp.multiplicative_replacement(skbio_comp.closure(x)))
    
    tmbed_df = pd.read_csv(tmbed_file, sep ='\t', comment='#')
    tmbed_df = tmbed_df[['prot_ID', 'P_B', 'P_H', 'P_S', 'P_i', 'P_o']]
    tmbed_df.set_index('prot_ID')
    tmbed_closed = tmbed_df.apply(closure_func, axis=1)
    tmbed_closed.reset_index()
    tmbed_closed['tmbed_slr'] = np.log((tmbed_closed['P_B'] + tmbed_closed['P_H']) / (tmbed_closed['P_S'] + tmbed_closed['P_i'] + tmbed_closed['P_o']))
    tmbed_closed = tmbed_closed[['prot_ID', 'tmbed_slr']]
    tmbed_closed.to_csv(out_file, sep='\t')
    return out_file

def process_tmbed_file(tmbed_file, ayu_outdir, keep_files = False):
    

    tmbed_comp_file = ayu_outdir +  'ayu.{}.tmbed_comp.tsv'.format(os.getpid())
    tmbed_slr_file = ayu_outdir +  'ayu.{}.tmbed_slr.tsv'.format(os.getpid())
    master_tmbed_record_list = parse_tmbed_file(tmbed_file)
    
    with open(tmbed_comp_file, 'w') as out_handle:
        result_header = None
        for finished_dict in get_avg_values(master_tmbed_record_list):
            if result_header is None:
                result_header = list(finished_dict.keys())
                out_handle.write('\t'.join(result_header) + '\n')
            out_handle.write('\t'.join([str(finished_dict[y]) for y in result_header]) + '\n')
    tmbed_slr_file = add_slr_to_tmbed(tmbed_comp_file, tmbed_slr_file)
    
    if not keep_files:
        os.remove(tmbed_comp_file)
    
    return tmbed_slr_file


#SignalP6
def parse_signalp6_file(infile, outfile):
    out_handle = open(outfile, 'w')
    out_handle.write('\t'.join(['prot_ID', 'sp_type', 'none_prob', 'sp_prob', 'cs_pos', 'cs_pos_prob']) + '\n')
    sp_type_list = ['NONE', 'SP', 'LIPO', 'TAT', 'TATLIPO', 'PILIN']
    with open(infile) as in_handle:
        for line in in_handle:
            if line[0] == '#': continue
            splitLine = line.rstrip('\n').split('\t')

            result_list = list(map(float, splitLine[2:-1]))
            max_index = result_list.index(max(result_list))
            max_value = splitLine[2 + max_index]
            sp_type = sp_type_list[max_index]

            #get max SP prob
            sp_result_list = list(map(float, splitLine[3:-1]))
            sp_max_index = sp_result_list.index(max(sp_result_list))
            sp_max_value = splitLine[3 + sp_max_index]

            if sp_type == 'NONE':
                out_handle.write('\t'.join([splitLine[0], sp_type, splitLine[2], sp_max_value, '0', '0']) + '\n')
            else:
                cs_pos_prob = splitLine[-1].split(':')[-1].strip(' ')
                cs_pos = splitLine[-1].split('.')[0].split(':')[-1].split('-')[0].strip(' ')
                out_handle.write('\t'.join([splitLine[0], sp_type, splitLine[2], sp_max_value, cs_pos, cs_pos_prob]) + '\n')
    out_handle.close()

    return outfile

def signalp6_logit_transform(in_file, out_file):
    def logit_func(x):
        x['sp_prob_logit'] = scipy.special.logit(x['sp_prob'])
        x['none_prob_logit'] = scipy.special.logit(x['none_prob'])
        return x

    signalp_df = pd.read_csv(in_file, sep='\t', comment='#')
    all_values = itertools.chain.from_iterable(signalp_df.drop(columns=['sp_type', 'cs_pos', 'cs_pos_prob', 'prot_ID']).values.tolist())
    all_values = sorted(list(set(all_values)))
    epsilon = (all_values[1]) / 2
    signalp_df['sp_prob'][signalp_df['sp_prob'] == 0] = epsilon
    signalp_df['sp_prob'][signalp_df['sp_prob'] >= 1] = 1-epsilon
    signalp_df = signalp_df.apply(logit_func, axis=1)
    signalp_df = signalp_df[['prot_ID', 'sp_prob_logit', 'none_prob_logit']]
    signalp_df.to_csv(out_file, mode = 'a', sep='\t')

    return out_file

def process_signalp_file(in_file, ayu_outdir, keep_files = False):
    signalp_proc_file = ayu_outdir + 'ayu.{}.signalp6.tsv'
    signalp_logit_file = ayu_outdir + 'ayu.{}.signalp6_logit.tsv'
    signalp_proc_file = parse_signalp6_file(in_file)
    signalp_logit_file = signalp6_logit_transform(signalp_proc_file, signalp_logit_file)
    if not keep_files:
        os.remove(signalp_proc_file)
    return signalp_logit_file


#IPC File
def parse_ipc2_file(in_file, out_file):
    out_handle = open(out_file, 'w')
    out_handle.write('prot_ID\tpI\n')
    with open(in_file) as in_handle:
        for id, seq in SimpleFastaParser(in_handle):
            new_id = id.split('||')[0]
            pi_value = float(id.split(':')[-1])
            out_handle.write('{}\t{}\n'.format(new_id, pi_value))
    out_handle.close()
    return out_file

def process_ipc2_file(in_file, ayu_outdir):
    ipc2_file = ayu_outdir + 'ayu.{}.ipc2.tsv'
    ipc2_file = parse_ipc2_file(in_file, ipc2_file)
    return ipc2_file


'''
def logit_func(x):
    return pd.Series(sp_special.logit(x))
def closure_func(x):
    return pd.Series(skbio_comp.multiplicative_replacement(skbio_comp.closure(x)))

    
    
tmbed_df.set_index('prot_ID', drop=True, inplace=True)
tmbed_colnames = tmbed_df.columns.values
tmbed_df2 = tmbed_df.apply(closure_func, axis=1)
tmbed_df2.rename({k:v for k,v in zip(tmbed_df2.columns.values, tmbed_colnames)}, axis='columns', inplace=True)
tmbed_df2.reset_index()                      
tmbed_df['tmbed_slr'] = np.log((tmbed_df2['P_B'] + tmbed_df2['P_H']) / (tmbed_df2['P_S'] + tmbed_df2['P_i'] + tmbed_df2['P_o']))

# https://stackoverflow.com/questions/43757977/replacing-values-greater-than-a-number-in-pandas-dataframe
# Apply this so we can do the logit transformation
# I should only apply logit transformation if I'm going to use this feature as an outcome. For this I don't give a shit, specially now
signalp_df.drop(columns=['sp_type', 'cs_pos', 'cs_pos_prob'], inplace=True)
signalp_df.set_index('prot_ID', inplace=True)
all_df_values = sorted(list(set(itertools.chain.from_iterable(signalp_df.values.tolist()))))
epsilon = (all_df_values[1]) / 2
signalp_df[signalp_df == 0] = epsilon
signalp_df[signalp_df >= 1] = 1-epsilon
signalp_df = signalp_df.apply(logit_func, axis=1)
'''




