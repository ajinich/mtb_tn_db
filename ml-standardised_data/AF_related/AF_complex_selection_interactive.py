import os
import pandas as pd
import seaborn as sns
import re
import matplotlib.pyplot as plt

colors_sns = sns.color_palette('colorblind')


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

            x = df_lfc[df_lfc.Rv_ID==x_rvid].values[0][2:]
            y = df_lfc[df_lfc.Rv_ID==y_rvid].values[0][2:]

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


# Load lfc data: 
fn_lfc = '../data/standardized_data/result_logfc_matrix_2021_11_15_BASIS_invitro.csv'
df_lfc = pd.read_csv(fn_lfc)
df_lfc.dropna(axis=0, inplace=True)
df_lfc.rename(columns={'rvid':'Rv_ID'}, inplace=True)


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



# Load GLS interaction data: 
fn = '/home/ajinich/Dropbox/KyuRhee/unknown_function/unknown_redox/data/GLS_TnSeq_v2/test_SI_data_1_fdr.001.xlsx'
df_interact = pd.read_excel(fn)
df_interact.sort_values(by = 'p_value_FDR', inplace=True)
print("number of unique pairs:", int(df_interact.shape[0]/2))


# Can we discard all the insertion sequences with strong artifacts?
str_insertion_1 = 'insertion sequence'
str_insertion_2 = 'transposase'
df_insertions = df_mtb_w_loc[df_mtb_w_loc['Protein names'].str.contains( '|'.join([str_insertion_1, str_insertion_2]))].copy()


list_rvid_insertions_1 = df_insertions.Rv_ID.tolist()
str_insertion_3 = df_mtb_w_loc[df_mtb_w_loc['Gene names'].str.contains('Rv2814')]['Gene names'].values[0]
list_rvid_insertions_2 = re.findall(re_str, str_insertion_3)

str_insertion_4 = df_mtb_w_loc[df_mtb_w_loc['Gene names'].str.contains('Rv3474')]['Gene names'].values[0]
list_rvid_insertions_3 = re.findall(re_str, str_insertion_4)

list_rvid_insertions = list(set(list_rvid_insertions_1 + list_rvid_insertions_2 + list_rvid_insertions_3))
df_interact = df_interact[ ~df_interact.lead_gene.isin(list_rvid_insertions) & ~df_interact.partner_gene.isin(list_rvid_insertions)] 

# Show scatter plots for top candidates: 
# Then you click and visually select which ones to prioritize: 
num_to_test = df_interact.shape[0]
df_test = df_interact.head(num_to_test).copy()


list_pairs = [sorted([x, y]) for x, y in zip(df_test['lead_gene'], df_test['partner_gene'])]
list_pairs_str = [g[0]+'_'+g[1] for g in list_pairs]
indexes = sorted([list_pairs_str.index(x) for x in set(list_pairs_str)])
list_pairs_unique = [list_pairs[i] for i in indexes]

indexes = sorted([list_pairs_str.index(x) for x in set(list_pairs_str)])
list_pairs_unique = [list_pairs[i] for i in indexes]
list_subset = []

list_keep = []
counter = 1
for list_rvid in list_pairs_unique:
	print(counter, 'out of', len(list_pairs_unique))
	len_width = int(len(list_rvid)*2.0)
	correlation_tile_plot(df_lfc, list_rvid, list_rvid, (len_width,len_width), colors_sns, dict_rvid_to_name, list_subset = list_subset, gene_names = False)
	# plt.tight_layout()
	plt.show()
	keep_or_toss = input("keep(y) or toss(n)?")
	if keep_or_toss == 'y':
		print('its a keeper!')
		list_keep.append(list_rvid)

	if counter%300 == 0:
		df_keep = pd.DataFrame()
		df_keep['Rv_ID_1'] = [rvs[0] for rvs in list_keep] 
		df_keep['Rv_ID_2'] = [rvs[1] for rvs in list_keep] 

		fn_out = '../data/AF_TnSeq_complex_keep_042822_'+str(counter)+'.xlsx'
		df_keep.to_excel(fn_out, index = False)
		list_keep = []

	
	counter += 1
	# plt.clf()

# # Store in a dataframe and write out to file. 


# # Rewrite so that you're storing temporarily == how would you do that? 
# df_keep = pd.DataFrame()
# df_keep['Rv_ID_1'] = [rvs[0] for rvs in list_keep] 
# df_keep['Rv_ID_2'] = [rvs[1] for rvs in list_keep] 

# fn_out = '../data/AF_TnSeq_complex_keep_042822.xlsx'
# df_keep.to_excel(fn_out, index = False)