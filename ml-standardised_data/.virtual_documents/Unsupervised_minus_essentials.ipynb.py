import pandas as pd
import os
import numpy as np


lfc = pd.read_csv('../data/standardized_data/result_logfc_matrix_2020_08_27.csv')


lfc.head()


lfc=lfc.dropna(axis=0)


mbio = pd.read_excel("../data/annotations/DeJesus_mbio.xlsx", header=1)
mbio.head()


print(lfc.shape)
lfc=pd.merge(lfc, mbio[['ORF ID', 'Name', 'Description', 'Final Call']], how='left', left_on='Rv_ID', right_on='ORF ID')
print(lfc.shape)


mcbwser = pd.read_excel("../data/annotations/Mycobacterium_tuberculosis_H37Rv_txt_v3.xlsx")
mcbwser.head()


mcbwser = mcbwser.drop_duplicates(subset=['Rv_ID'])


print(lfc.shape)
lfc=pd.merge(lfc, mcbwser[['Rv_ID', 'Functional_Category']], how='left', on='Rv_ID')
print(lfc.shape)


lfc=lfc.drop(columns='ORF ID')


lfc.head()


desc_columns=['Rv_ID', 'Name', 'Description', 'Final Call', 'Functional_Category']
value_cols=[col for col in lfc.columns if col not in desc_columns]
lfc=lfc[desc_columns + value_cols]


lfc.head()


#removing essentials
lfc = lfc[lfc['Final Call']get_ipython().getoutput("='ES']")




lfc.head()


from sklearn.decomposition import PCA
from plotnine import *


pca=PCA(3)
pca_results=pca.fit_transform(lfc[value_cols])
explained_var = pca.explained_variance_
pca_df=pd.DataFrame(pca_results, columns=['pca1', 'pca2', 'pca3'])


lfc=lfc.merge(pca_df, left_index=True, right_index=True)


lfc.head()


(ggplot(lfc, aes(x='pca1', y='pca2', color='Final Call')) 
 + geom_point() 
 + theme_light() 
 + xlab(f'pca1 {np.round(explained_var[0], 2)}') 
 + ylab(f'pca2 {np.round(explained_var[1], 2)}')
)


(ggplot(lfc, aes(x='pca1', y='pca2', color='Functional_Category'))
 +geom_point()
 +theme_light()
 +xlab(f'pca1 {np.round(explained_var[0], 2)}')
 +ylab(f'pca2 {np.round(explained_var[1], 2)}'))


from sklearn.manifold import TSNE


tsne = TSNE(n_components=2).fit_transform(lfc[value_cols])


tsne = pd.DataFrame(tsne, columns=['tsne1', 'tsne2'])


lfc=lfc.merge(tsne, left_index=True, right_index=True)


(ggplot(lfc, aes(x='tsne1', y='tsne2', color='Final Call'))+geom_point()+theme_light())


(ggplot(lfc, aes(x='tsne1', y='tsne2', color='Functional_Category'))+geom_point()+theme_light())


import scipy.cluster.hierarchy as shc
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 7))  
plt.title("Dendrograms")  
dend = shc.dendrogram(shc.linkage(lfc[value_cols], method='ward'))


from sklearn.cluster import AgglomerativeClustering 
hc = AgglomerativeClustering(n_clusters = 5, affinity = 'euclidean', linkage ='ward')
y_hc=hc.fit_predict(lfc[value_cols])


lfc['y_hc']=y_hc


lfc


(ggplot(lfc, aes(x='factor(y_hc)'))+geom_bar(aes(fill='Final Call'), position='fill')+theme_light()+ylab('Fractional Count'))


(ggplot(lfc, aes(x='factor(y_hc)'))+geom_bar(aes(fill='Functional_Category'), position='fill')+theme_light()+ylab('Fractional Count'))


lfc.to_csv('results/clustering_lfc_genes_minus_essentials.csv', index=False)



