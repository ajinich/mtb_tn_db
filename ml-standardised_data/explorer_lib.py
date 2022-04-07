import pandas as pd
import holoviews as hv 
import bokeh.io
import bokeh.plotting
import colorcet as cc
import os
import re
from sklearn.preprocessing import normalize
import matplotlib.pyplot as plt 
import torch
from sklearn.metrics import pairwise
import numpy as np


def import_interact_lfc_uniprot(norm):
    
    # interaction file
    path = '/home/ajinich/Dropbox/KyuRhee/unknown_function/unknown_redox/data/GLS_TnSeq_v2/'
    fn = 'test_SI_data_1_fdr.001.xlsx'
    fn_path = os.path.join(path, fn)
    df_interact = pd.read_excel(fn_path)

    # LFC dataset
    fn_lfc_basis = '../data/standardized_data/result_logfc_matrix_2021_11_15_BASIS_invitro.csv'
    df_lfc_basis = pd.read_csv(fn_lfc_basis)
    df_lfc_basis.dropna(axis=0, inplace=True)

    cols_data = df_lfc_basis.columns[1:]

    if norm: 
        X = df_lfc_basis[cols_data].values
        X_norm = normalize(X, norm='l2', axis=0)

        df_lfc_basis_norm_invitro = df_lfc_basis.copy()
        df_lfc_basis_norm_invitro[cols_data] = X_norm

        df_lfc = df_lfc_basis_norm_invitro.copy()
        
    else:
        df_lfc = df_lfc_basis.copy()

    # uniprot and dict_rvid_to_name
    fn = '/home/ajinich/Documents/repos/mtb_tn_db/data/annotations/uniprot_mtb_with_location.xlsx'
    df_mtb_w_loc = pd.read_excel(fn)
    df_mtb_w_loc = df_mtb_w_loc.fillna('')

    re_str = 'Rv\d\d\d\dc?'
    list_rvids = [re.findall(re_str, str_temp)[0] for str_temp in df_mtb_w_loc['Gene names']]
    df_mtb_w_loc['Rv_ID'] = list_rvids

    list_gene_names = [gn.split()[0] for gn in df_mtb_w_loc["Gene names"]]
    df_mtb_w_loc['gene_names'] = list_gene_names

    df_rvid_to_name = df_mtb_w_loc[['Rv_ID', 'gene_names']].copy() 

    dict_rvid_to_name = {}
    for index, row in df_rvid_to_name.iterrows():
        dict_rvid_to_name[row.Rv_ID] = row.gene_names


    return df_interact, df_lfc, df_mtb_w_loc, dict_rvid_to_name


def get_NN12(rvid_query, df_interact):
    # first nearest neighbors: 
    df_NN1 = df_interact[(df_interact.lead_gene==rvid_query) | (df_interact.partner_gene==rvid_query)].copy()
    list_rvid_NN1 = list(set(df_NN1.lead_gene.tolist() + df_NN1.partner_gene.tolist()))
    list_rvid_NN1.sort()

    # second nearest neighbors: 
    df_NN2 = df_interact[ (df_interact.lead_gene.isin(list_rvid_NN1)) | (df_interact.partner_gene.isin(list_rvid_NN1))].copy()
    list_rvid_NN2 = list(set(df_NN2.lead_gene.tolist() + df_NN2.partner_gene.tolist()))
    list_rvid_NN2.sort()
    
    # third nearest neighbors: 
    df_NN3 = df_interact[ (df_interact.lead_gene.isin(list_rvid_NN2)) | (df_interact.partner_gene.isin(list_rvid_NN2))].copy()
    list_rvid_NN3 = list(set(df_NN3.lead_gene.tolist() + df_NN3.partner_gene.tolist()))
    list_rvid_NN3.sort()
    
    # df_NN = pd.concat([df_NN1, df_NN2, df_NN3])
    # df_NN.drop_duplicates(inplace = True)

    return list_rvid_NN1, list_rvid_NN2 , list_rvid_NN3 #, df_NN


def correlation_tile_plot(df_lfc, list_rvid_x, list_rvid_y, fig_size, cols, dict_rvid_to_name, list_subset=[], gene_names = True ):
    
    if gene_names:
        list_gene_names_x = [dict_rvid_to_name[rvid] for rvid in list_rvid_x]
        list_gene_names_y = [dict_rvid_to_name[rvid] for rvid in list_rvid_y]
    else:
        list_gene_names_x = list_rvid_x
        list_gene_names_y = list_rvid_y
    
    fig, axs = plt.subplots(len(list_rvid_x), len(list_rvid_y), figsize=fig_size)
    if max([len(list_rvid_x), len(list_rvid_y)]) >= 10:
        FontSize = 20
        size_param = 40
    elif max([len(list_rvid_x), len(list_rvid_y)]) <= 5:
        FontSize = 20
        size_param = 70
    else:
        FontSize = 20
        size_param = 70
    for i in range(len(list_rvid_x)):
        for j in range(len(list_rvid_y)):
            x_rvid = list_rvid_x[i]
            y_rvid = list_rvid_y[j] 

            x = df_lfc[df_lfc.Rv_ID==x_rvid].values[0][1:]
            y = df_lfc[df_lfc.Rv_ID==y_rvid].values[0][1:]

            axs[i,j].scatter(x, y, s = size_param, alpha = 0.75, edgecolors='k', linewidths=3)
            axs[i,j].set_xticks([])
            axs[i,j].set_yticks([])
            if i == 0:
                axs[i,j].set_title(list_gene_names_y[j], fontsize = FontSize)
            
            if j == 0:
                axs[i,j].set_ylabel(list_gene_names_x[i], fontsize = FontSize)
            
            if list_rvid_x == list_rvid_y:
                if x_rvid in list_subset or y_rvid in list_subset:
                    axs[i,j].set_facecolor('xkcd:lightblue')
                if x_rvid in list_subset and y_rvid in list_subset:
                    axs[i,j].set_facecolor(cols[-2])
                if i==j:
                    axs[i,j].set_facecolor(cols[-3])

            # # testing grids for significant interactions: 
            # n_match = df_NN[ ((df_NN.lead_gene==x_rvid) & (df_NN.partner_gene==y_rvid)) | ((df_NN.lead_gene==y_rvid) & (df_NN.partner_gene==x_rvid))  ].shape[0]
            # if n_match:
            #     axs[i,j].grid(True)

def interactive_scatter_grid(df_lfc, list_rvid):

    cols_data = df_lfc.columns[1:].tolist()
    
    df_xy = df_lfc[ df_lfc.Rv_ID.isin(list_rvid) ][ ['Rv_ID']+cols_data ].copy()
    df_xy = df_xy.set_index('Rv_ID').T.rename_axis('screen').reset_index()

    list_hv = []
    for i in range(len(list_rvid)):
        hv_temp = df_xy.hvplot.scatter(x = list_rvid[i], y = list_rvid, width = 300, height = 300, size = 200, line_color='k', 
                                       line_width=3, hover_cols = ['screen'], subplots=True, fontsize = {'xlabel': '15pt'}, xlabel = list_rvid[i] ).cols(len(list_rvid))
        list_hv.append(hv_temp)

    return list_hv



def load_ESM_embeddings_and_df():
    # path to Mtb/M.smeg proteome: 
    fn = '/home/ajinich/Dropbox/KyuRhee/unknown_function/unknown_redox/data/mohammed/df_mtb_smeg_umap.csv'
    df_mtb_smeg = pd.read_csv(fn)
    df_mtb_smeg = df_mtb_smeg.fillna('')

    list_rvids = []
    re_str = 'Rv\d\d\d\dc?'
    for str_temp in df_mtb_smeg['Gene names']:
        re_match = re.findall(re_str, str_temp)
        if len(re_match):
            list_rvids.append(re_match[0])
        else:
            list_rvids.append('')
    df_mtb_smeg['Rv_ID'] = list_rvids

    # Path to ESM embeddings / representations: 
    path_rep = '/home/ajinich/Dropbox/KyuRhee/unknown_function/unknown_redox/data/mohammed/up_mtb_smeg_reprs/'
    EMB_LAYER = 33

    Xs_list = []
    list_err = []
    list_entries = df_mtb_smeg.Entry.tolist()

    for entry in list_entries:
        fn_full = os.path.join(path_rep, entry+'.pt')
        try:
            embs = torch.load(fn_full)
            Xs_list.append(embs['mean_representations'][EMB_LAYER])
        except:
            list_err.append(entry)
    X = torch.stack(Xs_list, dim=0).numpy()

    return X, df_mtb_smeg, list_entries 


def load_ESM_embeddings_and_df_10_prots():

    path_df = '/home/ajinich/Dropbox/KyuRhee/unknown_function/unknown_redox/data/mohammed/sakila_ESM/Proteomes/'
    list_fn = [os.path.join(path_df, fn) for fn in os.listdir(path_df)]
    list_df = [pd.read_csv(fn, sep='\t') for fn in list_fn]
    df_orgs = pd.concat(list_df, axis = 0)

    col = 'Gene names'
    df_orgs[col] = df_orgs[col].fillna('')
    list_rvids = []
    re_str = 'Rv\d\d\d\dc?'
    for str_temp in df_orgs[col]:
        re_match = re.findall(re_str, str_temp)
        if len(re_match):
            list_rvids.append(re_match[0])
        else:
            list_rvids.append('')
    df_orgs['Rv_ID'] = list_rvids

    path_rep = '/home/ajinich/Dropbox/KyuRhee/unknown_function/unknown_redox/data/mohammed/sakila_ESM/TotalProteomes/'
    EMB_LAYER = 33

    list_entries = [fn.split('.')[0] for fn in os.listdir(path_rep)]
    df_orgs = df_orgs[df_orgs.Entry.isin(list_entries)]
    df_orgs.reset_index(inplace=True, drop = True)

    Xs_list = []
    list_err = []
    list_entries = df_orgs.Entry.tolist()

    for entry in list_entries:
        fn_full = os.path.join(path_rep, entry+'.pt')
        try:
            embs = torch.load(fn_full)
            Xs_list.append(embs['mean_representations'][EMB_LAYER])
        except:
            list_err.append(entry)
    X = torch.stack(Xs_list, dim=0).numpy()

    return X, df_orgs, list_entries


def get_similar_prots(X, df, ind_query, list_entries, perc_th):

    x_test = X[ind_query,:]
    x_test = x_test.reshape(1, -1)

    mat_sim = pairwise.cosine_similarity(x_test, X)
    mat_d = 1 - mat_sim
    d_th = np.percentile(mat_d, perc_th)
    ind_th = list(np.where(mat_d <= d_th)[1])

    list_entries_th = [list_entries[i] for i in ind_th]

    df_th = pd.DataFrame()
    df_th['Entry'] = list_entries_th
    df_th['d']= mat_d[0, ind_th]
    df_th.sort_values(by = 'd', inplace=True)

    df_th = df_th.merge(df, on = 'Entry', how = 'left')

    return df_th