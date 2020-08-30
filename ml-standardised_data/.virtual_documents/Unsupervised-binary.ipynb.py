import pandas as pd
import os
import numpy as np


binary = pd.read_csv(
    '../data/standardized_data/result_bin_matrix_2020_08_27.csv')


binary.head()


binary = binary.dropna(axis=0)


mbio = pd.read_excel("../data/annotations/DeJesus_mbio.xlsx", header=1)
mbio.head()


print(binary.shape)
binary = pd.merge(binary, mbio[['ORF ID', 'Name', 'Description',
                                'Final Call']], how='left', left_on='Rv_ID', right_on='ORF ID')
print(binary.shape)


mcbwser = pd.read_excel(
    "../data/annotations/Mycobacterium_tuberculosis_H37Rv_txt_v3.xlsx")
mcbwser.head()


mcbwser = mcbwser.drop_duplicates(subset=['Rv_ID'])


print(binary.shape)
binary = pd.merge(
    binary, mcbwser[['Rv_ID', 'Functional_Category']], how='left', on='Rv_ID')
print(binary.shape)


binary = binary.drop(columns='ORF ID')





desc_columns = ['Rv_ID', 'Name', 'Description',
                'Final Call', 'Functional_Category']
value_cols = [col for col in binary.columns if col not in desc_columns]
binary = binary[desc_columns + value_cols]


binary[value_cols].sum(axis=1).value_counts()


binary['all_zero'] = binary.apply(
    lambda row: 'True' if row[value_cols].sum() == 0 else 'False', axis=1)


binary.head()


from sklearn.decomposition import PCA
from plotnine import *


pca = PCA(3)
pca_results = pca.fit_transform(binary[value_cols])
explained_var = pca.explained_variance_
pca_df = pd.DataFrame(pca_results, columns=['pca1', 'pca2', 'pca3'])


binary = binary.merge(pca_df, left_index=True, right_index=True)


(ggplot(binary, aes(x='pca1', y='pca2', color='Final Call'))
 + geom_point()
 + theme_light()
 + xlab(f'pca1 {np.round(explained_var[0], 2)}')
 + ylab(f'pca2 {np.round(explained_var[1], 2)}')
 )


(ggplot(binary, aes(x='pca1', y='pca2', color='Functional_Category'))
 + geom_point()
 + theme_light()
 + xlab(f'pca1 {np.round(explained_var[0], 2)}')
 + ylab(f'pca2 {np.round(explained_var[1], 2)}'))


(ggplot(binary, aes(x='pca1', y='pca2', color='all_zero'))
 + geom_point()
 + theme_light()
 + xlab(f'pca1 {np.round(explained_var[0], 2)}')
 + ylab(f'pca2 {np.round(explained_var[1], 2)}'))


from sklearn.manifold import TSNE


tsne = TSNE(n_components=2).fit_transform(binary[value_cols])


tsne = pd.DataFrame(tsne, columns=['tsne1', 'tsne2'])


binary = binary.merge(tsne, left_index=True, right_index=True)


(ggplot(binary, aes(x='tsne1', y='tsne2', color='Final Call'))+geom_point()+theme_light())


(ggplot(binary, aes(x='tsne1', y='tsne2', color='Functional_Category')) +
 geom_point()+theme_light())


(ggplot(binary, aes(x='tsne1', y='tsne2', color='all_zero'))+geom_point()+theme_light())


binary_minus_all_zeros = binary[binary.all_zero == 'False']
binary_minus_all_zeros = binary_minus_all_zeros[value_cols + desc_columns]
binary_minus_all_zeros.head()


tsne = TSNE(n_components=2).fit_transform(binary_minus_all_zeros[value_cols])


tsne = pd.DataFrame(tsne, columns=['tsne1', 'tsne2'])


binary_minus_all_zeros = binary_minus_all_zeros.merge(tsne, left_index=True, right_index=True)


(ggplot(binary_minus_all_zeros, aes(x='tsne1', y='tsne2', color='Final Call'))+geom_point()+theme_light())


(ggplot(binary_minus_all_zeros, aes(x='tsne1', y='tsne2', color='Functional_Category')) +
 geom_point()+theme_light())





import scipy.cluster.hierarchy as shc
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 7))
plt.title("Dendrograms")
dend = shc.dendrogram(shc.linkage(binary[value_cols], method='ward'))


from sklearn.cluster import AgglomerativeClustering
hc = AgglomerativeClustering(
    n_clusters=4, affinity='euclidean', linkage='ward')
y_hc = hc.fit_predict(binary[value_cols])


binary['y_hc'] = y_hc


binary


(ggplot(binary, aes(x='factor(y_hc)'))+geom_bar(aes(fill='Final Call'),
                                                position='fill')+theme_light()+ylab('Fractional Count'))


(ggplot(binary, aes(x='factor(y_hc)'))+geom_bar(aes(fill='Functional_Category'),
                                                position='fill')+theme_light()+ylab('Fractional Count'))


binary.to_csv('results/clustering_binary_genes.csv', index=False)
