import os
import pathlib

import numpy as np
import pandas as pd


lfc = pd.read_csv(
    '../data/standardized_data/result_logfc_matrix_2020_08_27.csv')
lfc.head()


qval = pd.read_csv(
    '../data/standardized_data/result_qval_matrix_2020_08_27.csv')
qval.head()


assert (lfc.columns == qval.columns).all()


lfc.isna().sum()


qval.isna().sum()


lfc[lfc['PE35_KO_vs_mbio_H37Rv'].isna() == True]['Rv_ID'].unique()


qval[qval['PE35_KO_vs_mbio_H37Rv'].isna() == True]['Rv_ID'].unique()


print(lfc.shape, qval.shape)
lfc = lfc.dropna(axis=0)
qval = qval.dropna(axis=0)
print(lfc.shape, qval.shape)


assert (lfc.columns == qval.columns).all()
assert (lfc.Rv_ID == qval.Rv_ID).all()


mcbwser = pd.read_excel(pathlib.Path.cwd().parents[0].joinpath(
    'data', 'annotations', 'Mycobacterium_tuberculosis_H37Rv_txt_v3.xlsx'))
mcbwser.head()


mcbwser[mcbwser['Rv_ID'].duplicated()]


mcbwser = mcbwser.drop_duplicates(subset=['Rv_ID'])


lfc_mb = pd.merge(lfc, mcbwser[['Rv_ID', 'Functional_Category']],
                  how='left', on='Rv_ID')
qval_mb = pd.merge(qval, mcbwser[['Rv_ID', 'Functional_Category']],
                   how='left', on='Rv_ID')


lfc_mb['Functional_Category'].value_counts(dropna=False)


lfc_mb.shape


qval_mb['Functional_Category'].value_counts(dropna=False)


qval_mb.shape


lfc_mb.to_csv('../data/standardized_data/cleaned_ML/lfc_mb.csv', index=False)
qval_mb.to_csv('../data/standardized_data/cleaned_ML/qval_mb.csv', index=False)


lfc_mb_filt = lfc_mb[~lfc_mb['Functional_Category'].isin(
    ['conserved hypotheticals', 'unknown'])]
qval_mb_filt = qval_mb[~qval_mb['Functional_Category'].isin(
    ['conserved hypotheticals', 'unknown'])]


assert (lfc_mb_filt.Rv_ID == qval_mb_filt.Rv_ID).all()


lfc_mb_filt.to_csv(
    '../data/standardized_data/cleaned_ML/lfc_mb_filt.csv', index=False)
qval_mb_filt.to_csv(
    '../data/standardized_data/cleaned_ML/qval_mb_filt.csv',  index=False)


data_cols = [col for col in qval_mb.columns if col not in [
    'Rv_ID', 'Functional_Category']]


bin_matrix_lfc = lfc_mb[data_cols].applymap(lambda x: x >= 1 or x <= -1)
bin_matrix_qval = qval_mb[data_cols].applymap(lambda x: x <= 0.05)


bin_matrix = bin_matrix_lfc & bin_matrix_qval
bin_matrix = bin_matrix.astype(int)


bin_matrix['Rv_ID'] = qval_mb['Rv_ID']
bin_matrix['Functional_Category'] = qval_mb['Functional_Category']


bin_matrix = bin_matrix[['Rv_ID'] +
                        [col for col in bin_matrix.columns if col get_ipython().getoutput("= 'Rv_ID']]")


bin_matrix.to_csv('../data/standardized_data/cleaned_ML/bin_mb.csv',  index=False)


bin_matrix_filt = bin_matrix[~bin_matrix['Functional_Category'].isin(
    ['conserved hypotheticals', 'unknown'])]


bin_matrix_filt.to_csv(
    '../data/standardized_data/cleaned_ML/bin_mb_filt.csv', index=False)



