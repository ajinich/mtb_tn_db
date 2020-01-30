import dash
import dash_core_components as dcc
import dash_daq as daq
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_table as dt
import pandas as pd
import numpy as np
from numpy import inf
import plotly.graph_objs as go
import os
import seaborn as sns
from random import random
import urllib.parse
import visdcc


main_data = pd.read_csv('Tn_library_DASH.csv', header=0)
# main_data = main_data.drop(columns='gene_name')
column_descriptors = pd.read_csv('column_descriptors.csv', header=0)
genes_data = pd.merge(main_data, column_descriptors, left_on='Expt', right_on='column_ID', how='left')
genes_data = genes_data.drop(columns=['sequence_data', 'column_ID', 'journal', 'first_author', 'last_author'])

cogs_df = pd.read_csv('all_cogs.csv', header=0)
cogs_desc = pd.read_csv('cog_names.csv', header=0, index_col=0, squeeze=True)
orf_details = pd.read_csv('ORF_details_final.csv')
unique_expts = main_data['Expt'].unique()
unique_Rvs = list(main_data['Rv_ID'].unique())
unique_genes= list(main_data['gene_name'].unique())
main_data['id'] = main_data['Rv_ID']
main_data.set_index('id', inplace=True, drop=False)

path_data = '../../data/'
df_uk = pd.read_csv(os.path.join(path_data, 'unknown_essentials/unknown_ALL_levels_essential_scores.csv'))
df_uk = df_uk[['Rv_ID', 'gene_name', 'UK_score_4']]





