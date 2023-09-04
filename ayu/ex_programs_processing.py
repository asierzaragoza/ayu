# This requires the prediction_results.txt file
import pandas as pd

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

def parse_tmbed_file(infile):
    """Parse TMBed file and store its info in a dict

    :param infile: TMBed file
    :type infile: str
    :return: List of dicts where key=prot_ID; values PRD: Predicted type, and the other values are the prediction values for each type.
    :rtype: list
    """ 
    tmbed_master_list = []

    with open(infile) as in_handle:
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
    
def get_tmbed_avg_values(tmbed_list, outfile):
    """_summary_

    :param tmbed_dict: TMBed list from parse_tmbed_file()
    :type tmbed_dict: dict
    :param outfile: file to store the average values.
    :type outfile: str
    """
    comment_line = '#B:Transmembrane beta strand; H:Transmembrane alpha helix; S:Signal peptide' + \
                   '; i:Non-transmembrane, inside; o:Non-transmembrane, outside\n'

    result_list = []
    for tmbed_dict in tmbed_list:
        tmbed_avg_dict = {'prot_ID':tmbed_dict['prot_ID']}
        for key in tmbed_dict:
            if key not in set(['prot_ID', 'PRD']):
                tmbed_avg_dict[key] = round(sum(tmbed_dict[key]) / len(tmbed_dict[key]), 4)
        result_list.append(tmbed_avg_dict) 
    
    with open(outfile) as out_handle:
        out_handle.write(comment_line)
        result_header = None
        for finished_dict in result_list:
            if result_header is None:
                result_header = list(finished_dict.keys())
                out_handle.write('\t'.join(result_header) + '\n')
            out_handle.write('\t'.join([str(finished_dict[y]) for y in result_header]) + '\n')
    
    return outfile


#todo: logit transformation to signalp values
#todo: SLR transformation to tmbed
#todo: add IPC2 info

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




