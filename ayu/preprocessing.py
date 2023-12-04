from Bio.SeqIO.FastaIO import SimpleFastaParser
import os
import pandas as pd
import pickle

# Get sequences
# give them an appropiate name, store old sequence names
# check for non-standard AA
# Check for length
# check for asterisk at the end!
# check for repeated sequence names

std_aa_set = set("ARNDCEQGHILKMFPSTWYV")

def remove_file(*files_to_remove):
    for filename in files_to_remove:
        try:
            os.remove(filename)
        except FileNotFoundError:
            pass


def process_fasta_files(file_name, out_file, rejected_file):
    out_handle = open(out_file, 'w')
    small_prots_set = set([])
    nonstd_prots_set = set([])
    duplicated_prots_list = []
    id_set = set([])
    with open(file_name) as in_handle:
        for id, seq in SimpleFastaParser(in_handle):
            # Remove stop codon
            if seq[-1] == '*':
                seq = seq[:-1]
            if len(seq) <= 50:
                small_prots_set.add(id)
                continue
            if len(std_aa_set.union(set(seq))) != len(std_aa_set):
                nonstd_prots_set.add(id)
                continue
            if len(id_set.union(set([id]))) == len(id_set):
                duplicated_prots_list.append(id)


            out_handle.write('>{}\n{}\n'.format(id, seq))
    out_handle.close()

    print(' finished!')
    if len(small_prots_set) + len(nonstd_prots_set) + len(duplicated_prots_list) > 0 :
        print('\t{} sequences are too small (<50 aa) and have been removed from further analyses.'.format(len(small_prots_set)))
        print('\t{} sequences contain nonstandard aa and have been removed from further analyses.'.format(len(nonstd_prots_set)))
        print('\t{} sequences have repeated headers and have been removed from further analyses.'.format(len(duplicated_prots_list)))
        with open(rejected_file, 'w') as out_handle:
            out_handle.write('Protein ID\tReason for rejection\n')
            for id in small_prots_set:
                out_handle.write('{}\tsmall size\n'.format(id))
            for id in nonstd_prots_set:
                out_handle.write('{}\tnon-standard aa\n'.format(id))
            for id in duplicated_prots_list:
                out_handle.write('{}\tduplicated header\n'.format(id))

    return out_file

def give_aliases(in_file, out_file, mapping_file):
    mapping_dict = {}
    out_handle = open(out_file, 'w')
    with open(in_file) as in_handle:
        seq_counter = 0
        pid = os.getpid()
        for id, seq in SimpleFastaParser(in_handle):
            new_id = 'seq_{}_{}'.format(pid, seq_counter)
            mapping_dict[new_id] = id
            out_handle.write('>{}\n{}\n'.format(new_id, seq))
            seq_counter += 1
    out_handle.close()
    
    with open(mapping_file, 'w') as out_handle:
        for key in mapping_dict:
            out_handle.write('{}\t{}\n'.format(key, mapping_dict[key]))
    return out_file, mapping_file

def divide_fasta_files(file_name, out_file_prefix, max_aa_size = 100000000):
    n_of_split_files = 0
    aa_added = 0
    divided_files_list = []
    out_filename = out_file_prefix + '__{}.faa'.format(n_of_split_files)
    out_handle = open(out_filename, 'w')
    with open(file_name) as in_handle:
        for id, seq in SimpleFastaParser(in_handle):
            out_handle.write('>{}\n{}\n'.format(id, seq))
            aa_added += len(seq)
            if aa_added >= max_aa_size:
                out_handle.close()
                divided_files_list.append(out_filename)
                n_of_split_files += 1
                aa_added = 0
                out_filename = out_file_prefix + '__{}.faa'.format(n_of_split_files)
                out_handle = open(out_filename, 'w')
    out_handle.close()
    divided_files_list.append(out_filename)
    n_of_split_files += 1
    
    return divided_files_list

def check_for_aliases(file_name, ayu_outdir, force_alias_change = True):
    #load alias file
    ayu_status_file, progress_status_dict = get_ayu_status(ayu_outdir)
    alias_file = progress_status_dict['alias_mapping']
    alias_dict = {}
    alias_set = set([])
    with open(alias_file) as in_handle:
        for line in in_handle:
            splitLine = line.rstrip('\n').split('\t')
            alias_dict[splitLine[1]] = splitLine[0]
            alias_set.add(splitLine[0])
    file_df = pd.read_csv(file_name, sep = '\t')
    file_id_set = set(file_df['prot_ID'].to_list())
    intersection_vs_original = file_id_set.intersection(set(alias_dict.keys()))
    intersection_vs_aliases = file_id_set.intersection(alias_set)
    #load rejected
    rejected_set = set([])
    with open(progress_status_dict['rejected_preproc']) as in_handle:
        in_handle.readline()
        for line in in_handle:
            splitLine = line.rstrip('\n').split('\t')
            rejected_set.add(splitLine[0])

    if len(intersection_vs_aliases) == len(file_id_set):
        print('No need for alias change')
        return file_name
    elif len(intersection_vs_original) == len(file_id_set):
        print('All headers found in alias file')
        file_df['prot_ID'] = file_df['prot_ID'].apply(lambda x: alias_dict[x])
        file_df.to_csv(file_name + '.AL', sep = '\t', index=None)
        remove_file(file_name)
        return file_name + '.AL'
    elif len(intersection_vs_original) > 0:
        non_found_seqs_wo_rejected = file_id_set - intersection_vs_original - rejected_set
        if len(non_found_seqs_wo_rejected) > 0:
            print('Some sequence IDs ({}) were not found in the alias file'.format(len(non_found_seqs_wo_rejected)))
        if force_alias_change:
            file_df = file_df[file_df['prot_ID'].isin(intersection_vs_original)]
            file_df['prot_ID'] = file_df['prot_ID'].apply(lambda x: alias_dict[x])
            file_df.to_csv(file_name + '.AL', sep = '\t', index=None)
            remove_file(file_name)
            return file_name + '.AL'
        else:
            return None
    else:
        print('No sequence IDs were not found in the alias file')
        return None
    
def find_ayu_status_file(ayu_outdir):
    ayu_file_list = [x for x in os.listdir(ayu_outdir) if x.split('.')[-1] == 'ayu']
    if len(ayu_file_list) == 1:
        return ayu_outdir + ayu_file_list[0]
    else:
        return None

def get_ayu_status(ayu_outdir):
    ayu_status_file = find_ayu_status_file(ayu_outdir)
    if ayu_status_file is None:
        print('Not able to find an ayu status file in folder {}. Please check the path and try again.')
        return (None, None)
    progress_status_dict = load_ayu_progress(ayu_status_file)
    return (ayu_status_file, progress_status_dict)

def load_ayu_progress(progress_file):
    progress_status_dict = None
    with open(progress_file, 'rb') as in_handle:
        progress_status_dict = pickle.load(in_handle)
    return progress_status_dict

def save_ayu_progress(progress_file, progress_status_dict):
    with open(progress_file, 'wb') as out_handle:
        pickle.dump(progress_status_dict, out_handle)

def get_seq_set(fasta_file):
    seq_set = set([])
    with open(fasta_file) as in_handle:
        for id, seq in SimpleFastaParser(in_handle):
            seq_set.add(id)
    return seq_set


def build_extr_df_for_merging(progress_status_dict, ayu_fasta_file, extr_type):
    #type = tmbed, signalp6, ipc2
    #get prot_ids from fasta_file to merge
    fasta_id_set = set([])

    ayu_file_prefix = '.'.join(ayu_fasta_file.split('.')[:-1])

    extr_file_dict = {'tmbed':'raw_tmbed_files', 'signalp6':'raw_signalp_files', 'ipc2':'raw_ipc2_files'}

    extr_file_list = progress_status_dict[extr_file_dict[extr_type]]
    with open(ayu_fasta_file) as in_handle:
        for id, seq in SimpleFastaParser(in_handle):
            fasta_id_set.add(id)
    
    extr_df_list = []
    for extr_file in extr_file_list:
        extr_df = pd.read_csv(extr_file, sep = '\t')
        extr_df = extr_df.loc[extr_df['prot_ID'].isin(fasta_id_set)]
        extr_df_list.append(extr_df)
    final_df = pd.concat(extr_df_list)

    del extr_df_list

    final_df.to_csv(ayu_file_prefix + '.{}.tsv'.format(extr_type), sep = '\t', index=False)
    return ayu_file_prefix + '.{}.tsv'.format(extr_type)


def merge_dfs(progress_status_dict):
    def merge_datasets(final_df, file_to_merge):
        dp_to_merge = pd.read_csv(file_to_merge, sep = '\t', comment = '#') 
        df_to_merge_size = len(dp_to_merge.index)
        final_df_size = len(final_df.index)
        if df_to_merge_size != final_df_size:
            print('Some protein IDs are not represented in file {}'.format(file_to_merge))
        final_df = final_df.merge(dp_to_merge, on='prot_ID', how='inner')
        return final_df
    result_list = []
    for fasta_file in progress_status_dict['final_ayu_files']:
        file_dict = progress_status_dict['final_ayu_files'][fasta_file]
        final_df = pd.read_csv(file_dict['ilr_aa_counts'], sep = '\t', comment='#')

        ayu_files_list = [
            file_dict['ilr_dp_counts'],
            file_dict['protein_costs'],
            file_dict['ilr_pqso'],
            file_dict['ilr_ppaac'],
            file_dict['transmemb'],
            file_dict['sp'],
            file_dict['pi']
        ]

        if None in ayu_files_list:
            list_names = ['DP counts', 'Protein costs', 'pQSO', 'pPAAC', 'TMBed', 'SignalP6', 'pI']
            indexes_w_none = [ayu_files_list.index(i) for i in list_names if i == None]
            print('The following features are not present: {}. Please load the corresponding files and try again.'.format(
                ';'.join([list_names[i] for i in indexes_w_none])))
            return None

        for ayu_file in ayu_files_list:
            final_df = merge_datasets(final_df, ayu_file)
        
        ayu_merged_file = '.'.join(fasta_file.split('.')[:-1]) + '.MRG'
        final_df.to_csv(ayu_merged_file, sep = '\t', index=False)
        result_list.append(ayu_merged_file)
    return result_list




