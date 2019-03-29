import pandas as pd
import os
import numpy as np

pd.options.mode.chained_assignment = None  # default='warn


def build_single_column_df(file_in, col_name, dir_name = '.'):
    # reads in SI excel file/table into single-column dataframe
    # works for SI tables that simply list genes called (conditionally) essential. 
    df_tn = pd.read_excel(os.path.join(dir_name, file_in))
    # this binarizes the data (column of 1's for all essentiality calls) 
    df_tn[col_name] = df_tn.shape[0]*[1]
    df_tn = df_tn[['Rv_ID', col_name]]
    return df_tn
    

def merge_with_whole_genome(df_WG, col_name, df_tn):
    # merges datafram with whole-genome Rv-ID + gene_name dataframe
    df_tn_WG = df_WG.merge(df_tn, how = 'left', on = 'Rv_ID')
    df_tn_WG[col_name].fillna(0, inplace = True)
   
    return df_tn_WG


def build_q_val_ratio_df(file_in, col_name, qval_col, log2FC_col, is_ratio_log2FC, dir_name = '.'):
    # reads in SI excel file/table into single-column dataframe
    df_tn = pd.read_excel(os.path.join(dir_name, file_in))
    
    target_qval_col = col_name+'_q_val'
    target_log2FC_col = col_name+'_log2FC'
    
    # q-values
    if qval_col != 'Nan': # if column exists: 
        df_tn[target_qval_col] = df_tn[qval_col]
    else:                       # else, just a bunch of NaN's
        df_tn[target_qval_col] = df_tn.shape[0]*[np.nan]
    
    # ratios / log2FC
    if is_ratio_log2FC == 'TRUE':
        df_tn[target_log2FC_col] = df_tn[log2FC_col]
    elif is_ratio_log2FC == 'FALSE':
        df_tn[target_log2FC_col] = np.log2(df_tn[log2FC_col])
    else:
        df_tn[target_log2FC_col] = df_tn.shape[0]*[np.nan]
    
    df_tn = df_tn[['Rv_ID', target_qval_col, target_log2FC_col]]
    
    return df_tn


def merge_with_whole_genome_qval_log2FC(df_WG, target_qval_col, target_log2FC_col, df_tn):
    # merges datafram with whole-genome Rv-ID + gene_name dataframe
    df_tn_WG = df_WG.merge(df_tn, how = 'left', on = 'Rv_ID')
    
    df_tn_WG[target_qval_col].fillna(np.nan, inplace = True)
    df_tn_WG[target_log2FC_col].fillna(np.nan, inplace = True)
   
    return df_tn_WG


def get_col_names(file, df_col_info):
    # what the q-value column is called in the original data
    qval_col = df_col_info[df_col_info.file == file].q_val_col_name.values[0]
    log2FC_col = df_col_info[df_col_info.file == file].ratio_col_name.values[0]
    is_ratio_log2FC = df_col_info[df_col_info.file == file].is_ratio_log2FC.values[0]

    # what we want to call it in the Tn-Matrix:
    col_name = df_col_info[df_col_info.file == file].col_name.values[0]
    target_qval_col = col_name+'_q_val'
    target_log2FC_col = col_name+'_log2FC'
    
    return qval_col, log2FC_col, is_ratio_log2FC, col_name, target_qval_col, target_log2FC_col


def get_qval_log2FC_func(file_list, df_col_info, df_WG, df_tn_qval_log2FC_ALL, dir_name = '.'):
    
    for file in file_list:
        print(file)
        # get source and target column names. 
        qval_col, log2FC_col, is_ratio_log2FC, col_name, target_qval_col, target_log2FC_col = get_col_names(file, df_col_info)
        # read the data.
        df_tn = pd.read_excel(os.path.join('../../data/Tn_datasets', file))
        
        # depending on whether there's a log2FC column or not: 
        if not(qval_col == 'NONE') and not(log2FC_col == 'NONE'):
            # get qval and logFC columns, and rename
            df_tn_qval_log2FC = df_tn[['Rv_ID', qval_col, log2FC_col]]
            df_tn_qval_log2FC.rename(columns = {qval_col:target_qval_col, log2FC_col:target_log2FC_col}, inplace=True)
            # deal with the format of the ratio column (is it already in log2FC format?)
            if not(is_ratio_log2FC):
                df_tn_qval_log2FC[target_log2FC_col] = np.log2(df_tn_qval_log2FC[target_log2FC_col])

        # here you're not adding a log2FC_col because there is no such data!
        elif not(qval_col == 'NONE') and log2FC_col == 'NONE':
            df_tn_qval_log2FC = df_tn[['Rv_ID', qval_col]]
            df_tn_qval_log2FC.rename(columns = {qval_col:target_qval_col}, inplace=True)
            # fill the log2FC column with NaN's
            df_tn_qval_log2FC[target_log2FC_col] = df_tn_qval_log2FC.shape[0]*[np.nan]

        # merge with full genome. 
        df_tn_qval_log2FC_WG = merge_with_whole_genome_qval_log2FC(df_WG, target_qval_col, target_log2FC_col, df_tn_qval_log2FC)

        # merge with full matrix of q-vals and log2FC changes:
        df_tn_qval_log2FC_ALL = df_tn_qval_log2FC_ALL.merge(df_tn_qval_log2FC_WG, how = 'inner', on = ['Rv_ID', 'gene_name'])

    return df_tn_qval_log2FC_ALL
