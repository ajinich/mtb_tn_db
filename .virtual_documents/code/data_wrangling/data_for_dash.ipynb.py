import pandas as pd
import numpy as np
import pathlib
import glob
import os


col_desc = pd.read_excel('../../data/column_descriptors_standardized.xlsx')


col_desc.head(10)


col_desc.column_ID_SI.unique()


# check if trackviews exist
# trackview_folders = list(pathlib.Path(
#     r'C:\Users\light\Dropbox\tuberculosis\tuberculosis\TrackView').glob('*'))
# trackview_names = [os.path.basename(x) for x in trackview_folders]


# all_conditions = list(col_desc['control']) + list(col_desc['experimental'])


# set(all_conditions).difference(set(trackview_names))


col_desc = pd.read_excel('../../data/column_descriptors_standardized.xlsx')


col_desc = col_desc.rename(columns = {'column_ID_2': 'column_ID_std'})


lfc = pd.read_csv(
    '../../data/standardized_data/result_logfc_matrix_2020_08_27.csv')


qval = pd.read_csv(
    '../../data/standardized_data/result_qval_matrix_2020_08_27.csv')


gene_names = pd.read_excel(
    '../../data/annotations/DeJesus_mbio.xlsx', header=1)
gene_names.head()


gene_names = gene_names.rename(
    columns={'Name': 'gene_name', 'ORF ID': 'Rv_ID'})


set(lfc.columns).difference(set(qval.columns))


set(qval.columns).difference(set(lfc.columns))


col_ID_standardized = list(col_desc.column_ID_std.unique())


set(qval.columns).difference(set(col_ID_standardized))


set(col_ID_standardized).difference(set(qval.columns))


lfc.head()


lfc = lfc.merge(gene_names[['Rv_ID', 'gene_name', 'Description']], how='left')


lfc = lfc.melt(id_vars=['Rv_ID', 'gene_name', 'Description'],
               var_name='Expt', value_name='log2FC')
qval = qval.melt(id_vars=['Rv_ID'], var_name='Expt', value_name='q-val')


print(lfc.shape, qval.shape)


lfc.head()


std_data = lfc.merge(qval, how='left', on=['Rv_ID', 'Expt'])


std_data.head()


std_data.isna().sum()


std_data['gene_name'] = std_data.apply(lambda x: x.Rv_ID if x.gene_name=='-' else x.gene_name, axis=1)
std_data['gene_name'] = std_data.apply(lambda x: x.Rv_ID if pd.isna(x.gene_name) else x.gene_name, axis=1)


std_data


std_data.isna().sum()


std_data = std_data[~std_data['log2FC'].isna()]


std_data.isna().sum()


std_data.to_csv('../../data/standardized_data_dash.csv', index=False)


# main_data.to_csv('../../data/main_data_dash.csv')


get_ipython().run_line_magic("reset", " -f")


import pandas as pd
import numpy as np
import pathlib
import glob
import os


col_desc = pd.read_excel('../../data/column_descriptors_standardized.xlsx')


col_desc = col_desc.rename(columns = {'column_ID_2': 'column_ID_std'})


lfc_si = pd.read_csv(
    '../../data/SI_datasets/SI_log2FC.csv')
lfc_si = lfc_si.drop(columns='gene_name')


qval_si = pd.read_csv(
    '../../data/SI_datasets/SI_qval.csv')
qval_si = qval_si.drop(columns='gene_name')


print(lfc_si.shape, qval_si.shape)


lfc_si


gene_names = pd.read_excel(
    '../../data/annotations/DeJesus_mbio.xlsx', header=1)
gene_names.head()


gene_names = gene_names.rename(
    columns={'Name': 'gene_name', 'ORF ID': 'Rv_ID'})


set(lfc_si.columns).difference(set(qval_si.columns))


set(qval_si.columns).difference(set(lfc_si.columns))


col_ID_SI = list(col_desc.column_ID_SI.unique())


set(qval_si.columns).difference(set(col_ID_SI))


set(col_ID_SI).difference(set(qval_si.columns))


lfc_si.head()


lfc_si = lfc_si.merge(gene_names[['Rv_ID', 'gene_name', 'Description']], how='left')


lfc_si.shape


lfc_si = lfc_si.melt(id_vars=['Rv_ID', 'gene_name', 'Description'],
               var_name='Expt', value_name='log2FC')
qval_si = qval_si.melt(id_vars=['Rv_ID'], var_name='Expt', value_name='q-val')


print(lfc_si.shape, qval_si.shape)


lfc_si.head()


si_data = lfc_si.merge(qval_si, how='left', on=['Rv_ID', 'Expt'])


si_data.head()


si_data.isna().sum()


si_data['gene_name'] = si_data.apply(lambda x: x.Rv_ID if x.gene_name=='-' else x.gene_name, axis=1)
si_data['gene_name'] = si_data.apply(lambda x: x.Rv_ID if pd.isna(x.gene_name) else x.gene_name, axis=1)


si_data


si_data.isna().sum()


# main_data.to_csv('../../data/main_data_dash.csv')


si_data.to_csv('../../data/si_data_dash.csv', index=False)


get_ipython().run_line_magic("reset", " -f")


import pandas as pd
import numpy as np
import pathlib
import glob
import os


col_desc = pd.read_excel('../../data/column_descriptors_standardized.xlsx')


col_desc = col_desc.rename(columns = {'column_ID_2': 'column_ID_std'})
col_desc = col_desc.drop(columns=['column_ID', 'wig_files'])


col_desc.head()


col_desc.head()


col_desc.to_csv('../../data/col_desc_dash.csv', index=False)
