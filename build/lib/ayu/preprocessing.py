from Bio.SeqIO.FastaIO import SimpleFastaParser
import os

# Get sequences
# give them an appropiate name, store old sequence names
# check for non-standard AA
# Check for length
# check for asterisk at the end!
# check for repeated sequence names

std_aa_set = set("ARNDCEQGHILKMFPSTWYV")

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
    return divided_files_list

def merge_dfs(df_list):
    pass