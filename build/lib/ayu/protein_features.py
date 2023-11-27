from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP
from importlib import resources as impresources
from . import data as ayu_data
import itertools
import functools
import json
import math
import numpy as np
import multiprocessing


#todo: Fix this! https://stackoverflow.com/questions/6028000/how-to-read-a-static-file-from-inside-a-python-package

data_file = (impresources.files(ayu_data) / 'protein_scales.json')
aa_alphabet = 'ACDEFGHIKLMNPQRSTVWY'
aaindex_dict = False
with open(data_file) as in_handle:
    aaindex_dict = json.loads(in_handle.read())

# AUXILIARY FUNCTIONS
def calculate_substitution(sequence, aascales_dict):
    value_list = []
    for letter in sequence:
        value_list.append(aascales_dict[letter])
    final_sub_value = sum(value_list) / len(value_list)
    return final_sub_value

def get_coup_factor_paac(seq, lag):
    """Calculate coupling factor as described in PMID: 11288174.
    :param seq: AA sequence
    :type seq: str
    :param lag: Correlation rank to calculate
    :type lag: int
    :return: Coupling factor
    :rtype: float
    """    
    def get_corr_values(aa1, aa2, scales_id):
        return math.pow(aaindex_dict[scales_id][aa2] - aaindex_dict[scales_id][aa1], 2)
    
    tau = 0
    for i in range((len(seq) - lag)):
        corr_factor_list = []
        aa1 = seq[i]
        aa2 = seq[i + lag]
        for scales_id in ['hydrophobicity_scales', 'hydrophilicity_scales', 'residuemass_scales']:
            corr_factor_list.append(get_corr_values(aa1, aa2, scales_id))
        tau += sum(corr_factor_list)/len(corr_factor_list)
    
    tau = tau / (len(seq) - lag)
    return round(tau, 4)

def get_socn_qso(prot_seq, lag):
    """Calculate coupling factor as described in PMID: 11097861.

    :param seq: AA sequence
    :type seq: str
    :param lag: Correlation rank to calculate
    :type lag: int
    :return: Coupling factor
    :rtype: float
    """    
    tau = 0
    for i in range((len(prot_seq) - lag)):
        aa1 = prot_seq[i]
        aa2 = prot_seq[i + lag]
        tau += math.pow(aaindex_dict['schneider_wrede_dp_dtce'][aa1 + aa2], 2)
    tau = tau / (len(prot_seq) - lag)
    return round(tau, 4)

def get_aa_counts(protseq):
    """Return amino acid counts.
    :param protseq: tuple(header, sequence), as returned by SimpleFastaParser.
    :type protseq: tuple
    :return: dict with amino acid counts and header (with key = 'prot_ID').
    :rtype: dict
    """    
    aa_comp = {k:0 for k in aa_alphabet}
    for item in protseq[1]:
        aa_comp[item] += 1
    
    new_aa_comp = {'AA_'+ k:v for k, v in aa_comp.items()}
    new_aa_comp['prot_ID'] = protseq[0]
    return new_aa_comp

def get_dp_counts(protseq):
    """Return dipeptide counts.
    :param protseq: tuple (header, sequence), as returned by SimpleFastaParser.
    :type protseq: tuple
    :return: dict with dipeptide counts and header (with key = 'prot_ID').
    :rtype: dict
    """    
    dp_comp = {''.join(k):0 for k in itertools.product(aa_alphabet, repeat=2)}
    for j in range(len(protseq[1]) - (len(protseq[1]) % 2)):
        current = protseq[1][j:j+2]
        if len(current) < 2: continue
        dp_comp[current] += 1
    new_dp_comp = {'DPept_'+ k:v for k, v in dp_comp.items()}
    new_dp_comp['prot_ID'] = protseq[0]
    return new_dp_comp

def get_protein_costs(protseq):
    """Calculate protein N/S/C average content, ATP cost and length.
    :param protseq: tuple(header, sequence), as returned by SimpleFastaParser.
    :type protseq: tuple
    :return: dict with protein costs and header (with key = 'prot_ID').
    :rtype: dict
    """    
    costs_dict = {'prot_ID':protseq[0]} 
    costs_dict['c_content'] = round(calculate_substitution(protseq[1], aaindex_dict['c_content']), 3)
    costs_dict['s_content'] = round(calculate_substitution(protseq[1], aaindex_dict['s_content']), 3)
    costs_dict['n_content'] = round(calculate_substitution(protseq[1], aaindex_dict['n_content']), 3)
    costs_dict['atp_cost'] = round(calculate_substitution(protseq[1], aaindex_dict['atp_cost']), 3)
    costs_dict['length'] = len(protseq[1])

    return costs_dict

def get_pi_biopython(protseq):
    """Calculate protein isoelectric point (pI) using Biopython. IPC2 is preferred.
    :param protseq: tuple(header, sequence), as returned by SimpleFastaParser.
    :type protseq: tuple
    :return: dict with protein costs and header (with key = 'prot_ID').
    :rtype: dict
    """   
    pi_dict = {'prot_ID':protseq[0]}
    pi_dict = {'pI':IP(protseq[1]).pi()}

    return pi_dict

def get_ppaac(prot_seq, max_lag):
    """Calculate partial Pseudo Amino Acid Composition (pPAAC).
    :param protseq: tuple (header, sequence), as returned by SimpleFastaParser.
    :type protseq: tuple
    :param max_lag: Number of rank correlations to calculate. Cannot be higher than protein length.
    :type max_lag: int
    :return: dict with  max_lag pPAAC values and header (with key = 'prot_ID').
    :rtype: dict
    """    
    result_dict = {'prot_ID':prot_seq[0]}
    assert(len(prot_seq[1]) > max_lag)
    tau_list = []
    for i in range(1, max_lag + 1):
        tau_list.append(get_coup_factor_paac(prot_seq[1], i))
    total_tau = sum(tau_list)
    #Sometimes this happens? I still don't know why
    if total_tau == 0:
        print('ERROR: protein ID {} has encountered an error in pPAAC calculations'.format(prot_seq[0]))
        return None
    for i in range(len(tau_list)):
        result_dict['pPseAAC_{}'.format(i+1)] = round(tau_list[i] / total_tau, 4)
    
    return result_dict

def get_pqso(prot_seq, max_lag):
    """Calculate partial Quasi Sequence Order (pQSO).
    :param protseq: tuple (header, sequence), as returned by SimpleFastaParser.
    :type protseq: tuple
    :param max_lag: Number of rank correlations to calculate. Cannot be higher than protein length.
    :type max_lag: int
    :return: dict with  max_lag pQSO values and header (with key = 'prot_ID').
    :rtype: dict
    """    
    result_dict = {'prot_ID':prot_seq[0]}
    assert(len(prot_seq[1]) > max_lag)
    tau_list = []
    for i in range(1, max_lag + 1):
        tau_list.append(get_socn_qso(prot_seq[1], i))
    total_tau = sum(tau_list)
    # Sometimes this happens? I still don't know why
    if total_tau == 0:
        print('ERROR: protein ID {} has encountered an error in pQSO calculations'.format(prot_seq[0]))
        return None
    for i in range(len(tau_list)):
        result_dict['pQSO_{}'.format(i+1)] = round(tau_list[i] / total_tau, 4)
    
    return result_dict

#MAIN FUNCTIONS
def process_aa_counts(input_file, output_file, n_threads):
    pool = multiprocessing.Pool(n_threads)
    in_handle = open(input_file)
    out_handle = open(output_file, 'w')
    result_header = None

    for result in pool.imap_unordered(get_aa_counts, SimpleFastaParser(in_handle)):
        if result is None: continue
        if result_header is None:
            result_header = sorted(list(result.keys()))
            out_handle.write('\t'.join(result_header) + '\n')
        out_handle.write('\t'.join([str(result[y]) for y in result_header]) + '\n')  
    out_handle.close()
    return output_file
        
def process_dp_counts(input_file, output_file, n_threads):
    pool = multiprocessing.Pool(n_threads)
    in_handle = open(input_file)
    out_handle = open(output_file, 'w')
    result_header = None

    for result in pool.imap_unordered(get_dp_counts, SimpleFastaParser(in_handle)):
        if result is None: continue
        if result_header is None:
            result_header = sorted(list(result.keys()))
            out_handle.write('\t'.join(result_header) + '\n')
        out_handle.write('\t'.join([str(result[y]) for y in result_header]) + '\n') 
    out_handle.close()
    return output_file

def process_protein_costs(input_file, output_file, n_threads):
    pool = multiprocessing.Pool(n_threads)
    in_handle = open(input_file)
    out_handle = open(output_file, 'w')
    result_header = None

    for result in pool.imap_unordered(get_protein_costs, SimpleFastaParser(in_handle)):
        if result is None: continue
        if result_header is None:
            result_header = sorted(list(result.keys()))
            out_handle.write('\t'.join(result_header) + '\n')
        out_handle.write('\t'.join([str(result[y]) for y in result_header]) + '\n') 
    out_handle.close()
    return output_file

def process_pqso(input_file, output_file, n_threads, max_lag):
    pool = multiprocessing.Pool(n_threads)
    in_handle = open(input_file)
    out_handle = open(output_file, 'w')
    result_header = None
    partial_func_qso = functools.partial(get_pqso, max_lag = max_lag)
    for result in pool.imap_unordered(partial_func_qso, SimpleFastaParser(in_handle)):
        if result is None: continue
        if result_header is None:
            result_header = sorted(list(result.keys()))
            out_handle.write('\t'.join(result_header) + '\n')
        out_handle.write('\t'.join([str(result[y]) for y in result_header]) + '\n')   
    out_handle.close()
    return output_file

def process_ppaac(input_file, output_file, n_threads, max_lag):
    pool = multiprocessing.Pool(n_threads)
    in_handle = open(input_file)
    out_handle = open(output_file, 'w')
    result_header = None
    partial_func_paac = functools.partial(get_ppaac, max_lag = max_lag)
    for result in pool.imap_unordered(partial_func_paac, SimpleFastaParser(in_handle)):
        if result is None: continue
        if result_header is None:
            result_header = sorted(list(result.keys()))
            out_handle.write('\t'.join(result_header) + '\n')
        out_handle.write('\t'.join([str(result[y]) for y in result_header]) + '\n')   
    out_handle.close()
    return output_file

def process_pi(input_file, output_file, n_threads):
    pool = multiprocessing.Pool(n_threads)
    in_handle = open(input_file)
    out_handle = open(output_file, 'w')
    result_header = None

    for result in pool.imap_unordered(get_pi_biopython, SimpleFastaParser(in_handle)):
        if result is None: continue
        if result_header is None:
            result_header = sorted(list(result.keys()))
            out_handle.write('\t'.join(result_header) + '\n')
        out_handle.write('\t'.join([str(result[y]) for y in result_header]) + '\n') 
    out_handle.close()
    return output_file