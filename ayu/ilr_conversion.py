import skbio.stats.composition as skbio_comp
import sklearn.decomposition as sk_decomp
import sklearn.preprocessing as sk_prepro
import scipy.special as sp_special
import numpy as np
import pandas as pd
import os
import multiprocessing

def closure_func(x):
    return pd.Series(skbio_comp.multiplicative_replacement(skbio_comp.closure(x)))

def ilr_func(x):
    return pd.Series(skbio_comp.ilr(x))

def closure_func_parallel(df):
    return df.apply(closure_func, axis = 1)

def ilr_func_parallel(df):
    return df.apply(closure_func, axis = 1)

def parallelize_dataframe_apply(df, n_threads, func):
    df_split = np.array_split(df, n_threads)
    with multiprocessing.Pool(n_threads) as p:
        df = pd.concat(p.map(func, df_split))
    return df

def process_file_parallel(filename, n_threads):
    out_comp_file = '.'.join(filename.split('.')[:-1]) + '.comp.tsv'
    out_ilr_file = '.'.join(filename.split('.')[:-1]) + '.ilr.tsv'
    #Get and store comment lines
    comment_lines = ''
    with open(filename) as in_handle:
        for line in in_handle:
            if line[0] == '#':
                comment_lines += line
    
    aa_comp_df = pd.read_csv(filename, sep = '\t', comment = '#')
    aa_comp_df.set_index('prot_ID', drop=True, inplace=True)
    
    closure_df = parallelize_dataframe_apply(aa_comp_df, n_threads, closure_func_parallel)
    closure_df.rename({k:v for k,v in zip(closure_df.columns.values, 
                        ['Comp_' + str(i) for i in aa_comp_df.columns.values]) }, axis='columns', inplace=True)
    
    comment_lines_fst = comment_lines + '#Applied Composition + Closure\n'
    with open(out_comp_file, 'w') as out_handle:
        out_handle.write(comment_lines_fst)
    closure_df.to_csv(out_comp_file, mode = 'a', sep='\t')

    closure_columns = None
    if len(closure_df.columns.values[0].split('_') <= 2):
        closure_columns = [x.split('_')[0] for x in closure_df.columns.values]
    else:
        closure_columns = ['_'.join(x.split('_')[:-1]) for x in closure_df.columns.values]
    
    ilr_aa_comp_df = parallelize_dataframe_apply(closure_df, n_threads, ilr_func_parallel)
    ilr_aa_comp_df.rename({k:v for k,v in zip(range(len(ilr_aa_comp_df.columns) + 1), 
                ['ILR_' + i + '_' + str(j) for i in closure_columns for j in range(len(ilr_aa_comp_df.columns) + 1)])}, axis='columns', inplace=True)
    
    comment_lines += '#Applied ILR transformation\n'
    with open(out_ilr_file, 'w') as out_handle:
        out_handle.write(comment_lines)
    ilr_aa_comp_df.to_csv(out_ilr_file, mode='a', sep='\t')
    return out_ilr_file

def process_file(filename):
    out_comp_file = '.'.join(filename.split('.')[:-1]) + '.comp.tsv'
    out_ilr_file = '.'.join(filename.split('.')[:-1]) + '.ilr.tsv'
    #Get and store comment lines
    comment_lines = ''
    with open(filename) as in_handle:
        for line in in_handle:
            if line[0] == '#':
                comment_lines += line
    
    aa_comp_df = pd.read_csv(filename, sep = '\t', comment = '#')
    aa_comp_df.set_index('prot_ID', drop=True, inplace=True)
    
    closure_df = aa_comp_df.apply(closure_func, axis=1)
    closure_df.rename({k:v for k,v in zip(closure_df.columns.values, 
                        ['Comp_' + str(i) for i in aa_comp_df.columns.values]) }, axis='columns', inplace=True)
    
    comment_lines_fst = comment_lines + '#Applied Composition + Closure\n'
    with open(out_comp_file, 'w') as out_handle:
        out_handle.write(comment_lines_fst)
    closure_df.to_csv(out_comp_file, mode = 'a', sep='\t')

    closure_columns = ['_'.join(x.plit('_')[1:-1]) for x in closure_df.columns.values]
    ilr_aa_comp_df = closure_df.apply(ilr_func, axis=1)
    ilr_aa_comp_df.rename({k:v for k,v in zip(range(len(ilr_aa_comp_df.columns) + 1), 
                ['ILR_' + i + '_' + str(j) for i in closure_columns for j in range(len(ilr_aa_comp_df.columns) + 1)])}, axis='columns', inplace=True)
    
    comment_lines += '#Applied Composition + Closure + ILR transformation\n'
    with open(out_ilr_file, 'w') as out_handle:
        out_handle.write(comment_lines)
    ilr_aa_comp_df.to_csv(out_ilr_file, mode='a', sep='\t')
    return filename