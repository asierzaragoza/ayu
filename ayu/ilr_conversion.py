import skbio.stats.composition as skbio_comp
import sklearn.decomposition as sk_decomp
import sklearn.preprocessing as sk_prepro
import scipy.special as sp_special
import numpy as np
import pandas as pd
import functools
import os
import multiprocessing

def closure_func(x):
    return pd.Series(skbio_comp.multiplicative_replacement(skbio_comp.closure(x)))

def ilr_func(x):
    return pd.Series(skbio_comp.ilr(x))

def closure_func_parallel(df):
    return df.apply(closure_func, axis = 1)

def ilr_func_parallel(df):
    return df.apply(ilr_func, axis = 1)

def parallelize_dataframe_apply(df, n_threads, func):
    df_split = np.array_split(df, n_threads)
    with multiprocessing.Pool(n_threads) as p:
        df = pd.concat(p.map(func, df_split))
    return df

def parallelize_dataframe_apply2(input_df, func):
    #input_df.set_index('prot_ID', drop=True, inplace=True)
    # Remove all rows that only have 0s
    input_df = input_df.loc[(input_df!=0).any(axis=1)]
    output_df = func(input_df)
    # Reduce the number of decimal numbers so the df is not huge
    output_df = output_df.round(decimals=6)
    return output_df



def process_file_parallel2_closure(filename, n_threads, chunksize = 250000):
    out_comp_file = '.'.join(filename.split('.')[:-1]) + '.comp.tsv'
    out_ilr_file = '.'.join(filename.split('.')[:-1]) + '.ilr.tsv'
    #Get and store comment lines
    comment_lines = ''
    header = None
    with open(filename) as in_handle:
        for line in in_handle:
            if line[0] == '#':
                comment_lines += line
    #Write comment lines
    comment_lines += '#Applied Composition + Closure\n'
    with open(out_comp_file, 'w') as out_handle:
        out_handle.write(comment_lines)
    
    input_columns = pd.read_csv(filename, sep = '\t', comment='#', nrows=10, index_col='prot_ID').columns.values
    input_df_chunks = pd.read_csv(filename, sep = '\t', comment='#', index_col='prot_ID', chunksize=chunksize)
    pool = multiprocessing.Pool(n_threads)
    closure_partial_func = functools.partial(parallelize_dataframe_apply2, func=closure_func_parallel)
    for result_df in pool.imap_unordered(closure_partial_func, input_df_chunks):
        if header is None:
            header = tuple(result_df.columns.values)
            result_df.rename({k:v for k,v in zip(result_df.columns.values, 
                        ['Comp_' + str(i) for i in input_columns]) }, axis='columns', inplace=True)
            result_df.to_csv(out_comp_file, mode = 'a', sep = '\t')
        
        else:
            result_df = result_df[header]
            result_df.to_csv(out_comp_file, mode = 'a', sep = '\t', header=False)
    
    return out_comp_file

def process_file_parallel2_ilr(filename, n_threads, chunksize = 250000):
    out_ilr_file = '.'.join(filename.split('.')[:-1]) + '.ilr.tsv'
    #get the header
    closure_columns = pd.read_csv(filename, sep = '\t', comment='#', nrows=10, index_col='prot_ID').columns.values
    
    if len(closure_columns[0].split('_')) <= 2:
        closure_columns = [x.split('_')[0] for x in closure_columns]
    else:
        closure_columns = ['_'.join(x.split('_')[:-1]) for x in closure_columns]
    
    #Get and store comment lines
    comment_lines = ''
    header = None
    with open(filename) as in_handle:
        for line in in_handle:
            if line[0] == '#':
                comment_lines += line
    #Write comment lines
    comment_lines += '#Applied ILR transformation\n'
    with open(out_ilr_file, 'w') as out_handle:
        out_handle.write(comment_lines)

    input_df_chunks = pd.read_csv(filename, sep = '\t', comment='#', index_col='prot_ID', chunksize=chunksize)
    pool = multiprocessing.Pool(n_threads)
    ilr_partial_func = functools.partial(parallelize_dataframe_apply2, func=ilr_func_parallel)
    for result_df in pool.imap_unordered(ilr_partial_func, input_df_chunks):
        if header is None:
            header = tuple(result_df.columns.values)
            result_df.rename({k:v for k,v in zip(range(len(result_df.columns) + 1), 
                ['ILR_' + i + '_' + str(j) for i in closure_columns for j in range(len(result_df.columns) + 1)])},
                axis='columns', inplace=True)
            result_df.to_csv(out_ilr_file, mode = 'a', sep = '\t')
        else:
            result_df = result_df[header]
            result_df.to_csv(out_ilr_file, mode = 'a', sep = '\t', header=False)
    
    return out_ilr_file



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
    if len(closure_df.columns.values[0].split('_')) <= 2:
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