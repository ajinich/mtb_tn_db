import pandas as pd
import numpy as np
import os


lfc = pd.read_csv('../data/standardized_data/result_logfc_matrix_2020_08_27.csv')


lfc.head()


## REMOVE ESSENTIALS!
mbio = pd.read_excel("../data/annotations/DeJesus_mbio.xlsx", header=1)
mbio.head()


print(lfc.shape)
lfc = pd.merge(lfc, mbio[['ORF ID', 'Final Call']], how='left', left_on='Rv_ID', right_on='ORF ID')
print(lfc.shape)


lfc = lfc[lfc['Final Call'] get_ipython().getoutput("= 'ES']")
lfc = lfc.drop(columns = ['ORF ID', 'Final Call'])
print(lfc.shape)


lfc=lfc.dropna(axis=0)
lfc=lfc.T.reset_index()
lfc.columns=list(lfc.iloc[0,:])
lfc=lfc.iloc[1:,:]
lfc=lfc.rename(columns={'Rv_ID': 'Dataset'})
lfc=lfc.set_index('Dataset')


lfc.head()


value_cols=[col for col in lfc.columns if colget_ipython().getoutput("='Dataset']")


from sklearn.decomposition import PCA
from plotnine import *


pca=PCA(3)
pca_results=pca.fit_transform(lfc[value_cols])
# print(pca_results)
explained_var = pca.explained_variance_ratio_
print(pca.explained_variance_ratio_)
pca_df=pd.DataFrame(pca_results, columns=['pca1', 'pca2', 'pca3'])
pca_df=pca_df.set_index(lfc.index)


lfc=lfc.merge(pca_df, left_index=True, right_index=True)


lfc


(ggplot(lfc, aes(x='pca1', y='pca2')) 
 + geom_point() 
 + theme_light() 
 + geom_label(label=lfc.index, size=5, nudge_x=10)
 + xlab(f'pca1 {np.round(explained_var[0], 2)}') 
 + ylab(f'pca2 {np.round(explained_var[1], 2)}')
)


from sklearn.manifold import TSNE


tsne = TSNE(n_components=2).fit_transform(lfc[value_cols])


tsne = pd.DataFrame(tsne, columns=['tsne1', 'tsne2'])
tsne.index=lfc.index


lfc=lfc.merge(tsne, left_index=True, right_index=True)


(ggplot(lfc, aes(x='tsne1', y='tsne2'))+geom_point()+theme_light())


lfc.columns.name = None


import scipy.cluster.hierarchy as shc
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(15, 14))
ax = fig.add_subplot(1, 1, 1)
shc.dendrogram(shc.linkage(lfc[value_cols], method='ward'),  color_threshold=150,labels=lfc.index)
ax.tick_params(axis='x', which='major', labelsize=13)
ax.tick_params(axis='y', which='major', labelsize=13)
fig.savefig('dataset_dendrogram.jpg', bbox_inches = 'tight', dpi=500)



from sklearn.cluster import AgglomerativeClustering 
hc = AgglomerativeClustering(n_clusters = 3, affinity = 'euclidean', linkage ='ward')
y_hc=hc.fit_predict(lfc[value_cols])


lfc['y_hc']=y_hc


lfc
