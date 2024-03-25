import pandas as pd
import numpy as np
import os


lfc = pd.read_csv(
    '../data/standardized_data/result_logfc_matrix_2020_08_27.csv')


lfc.head()


# REMOVE ESSENTIALS!
# mbio = pd.read_excel("../data/annotations/DeJesus_mbio.xlsx", header=1)
# mbio.head()


# print(lfc.shape)
# lfc = pd.merge(lfc, mbio[['ORF ID', 'Final Call']], how='left', left_on='Rv_ID', right_on='ORF ID')
# print(lfc.shape)


# lfc = lfc[lfc['Final Call'] get_ipython().getoutput("= 'ES']")
# lfc = lfc.drop(columns = ['ORF ID', 'Final Call'])
# print(lfc.shape)


lfc.isna().sum()


lfc = lfc.dropna(axis=0)
lfc = lfc.T.reset_index()
lfc.columns = list(lfc.iloc[0, :])
lfc = lfc.iloc[1:, :]
lfc = lfc.rename(columns={'Rv_ID': 'Dataset'})


dataset_labels = pd.read_csv('../data/standardized_data_labels.csv')
lfc = lfc.merge(dataset_labels, how='left', on='Dataset')
lfc.head()


value_cols = [col for col in lfc.columns if col not in ['Dataset', 'Label']]


from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from plotnine import *


lfc_scaled = StandardScaler().fit_transform(lfc[value_cols])
pca = PCA(3)
pca_results = pca.fit_transform(lfc_scaled)
# print(pca_results)
explained_var = pca.explained_variance_ratio_
print(pca.explained_variance_ratio_)
pca_df = pd.DataFrame(pca_results, columns=['pca1', 'pca2', 'pca3'])


lfc = lfc.merge(pca_df, left_index=True, right_index=True)


lfc


(ggplot(lfc, aes(x='pca1', y='pca2', color='Label'))
 + geom_point()
 + theme_light()
 #+ geom_label(label=lfc.index, size=5, nudge_x=10)
 + xlab(f'pca1 {np.round(explained_var[0], 2)}')
 + ylab(f'pca2 {np.round(explained_var[1], 2)}')
 )


(ggplot(lfc, aes(x='pca2', y='pca3', color='Label'))
 + geom_point()
 + theme_light()
 #+ geom_label(label=lfc.index, size=5, nudge_x=10)
 + xlab(f'pca2 {np.round(explained_var[1], 2)}')
 + ylab(f'pca3 {np.round(explained_var[2], 2)}')
 )


from sklearn.manifold import TSNE


tsne = TSNE(n_components=2).fit_transform(lfc_scaled)


tsne = pd.DataFrame(tsne, columns=['tsne1', 'tsne2'])


lfc = lfc.merge(tsne, left_index=True, right_index=True)


(ggplot(lfc, aes(x='tsne1', y='tsne2'))+geom_point()+theme_light())


import scipy.cluster.hierarchy as shc
import palettable
# shc.set_link_color_palette(['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf'])
import matplotlib.pyplot as plt
Z = shc.linkage(lfc[value_cols], method='ward', optimal_ordering=True)


# shc.set_link_color_palette(['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666'])
# shc.set_link_color_palette(['#000000']*9)
palette = palettable.colorbrewer.qualitative.Dark2_8.hex_colors
palette.append('#000000')
shc.set_link_color_palette(palette)


lfc.Label.isna().sum()


print(lfc[value_cols].shape, Z.shape)





fig = plt.figure(figsize=(8, 14))
ax = fig.add_subplot(1, 1, 1)
shc.dendrogram(Z,  color_threshold=150, labels=lfc.Dataset.to_list(),
               orientation='left', above_threshold_color='grey')
ax.tick_params(axis='x', which='major', labelsize=11)
ax.tick_params(axis='y', which='major', labelsize=10.5)
plt.axvline(x=150, c='grey', lw=2, linestyle='dashed')
# label_palette = palettable.matplotlib.Inferno_13.hex_colors

# # transforme the 'cyl' column in a categorical variable. It will allow to put one color on each level.
# cat_label = pd.Categorical(lfc['Label'])
# # print(cat_label.codes)
# my_color = cat_label.codes

# # Apply the right color to each label
# ax = plt.gca()
# xlbls = ax.get_ymajorticklabels()
# num = -1
# for lbl in xlbls:
#     num += 1
#     val = my_color[num]
#     lbl.set_color(label_palette[val])


fig.savefig('dataset_dendrogram.jpg', bbox_inches='tight', dpi=500)





# from sklearn.cluster import AgglomerativeClustering
#  hc = AgglomerativeClustering(n_clusters = 3, affinity = 'euclidean', linkage ='ward')
# y_hc=hc.fit_predict(lfc[value_cols])


# lfc['y_hc']=y_hc


lfc


from sklearn.cluster import KMeans


k = lfc.Label.nunique()


kmeans = KMeans(n_clusters=k, random_state=42).fit(lfc[value_cols])


y_pred = kmeans.predict(lfc[value_cols])


y_pred


lfc.Label


from sklearn.metrics import adjusted_rand_score


adjusted_rand_score(y_pred, lfc.Label)


from sklearn.cluster import AgglomerativeClustering
hc = AgglomerativeClustering(
    n_clusters=k, affinity='euclidean', linkage='ward')
y_hc = hc.fit_predict(lfc[value_cols])


adjusted_rand_score(y_hc, lfc.Label)
