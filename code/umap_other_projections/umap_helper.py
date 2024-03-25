
import numpy as np
import tqdm
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import scipy.stats as st
import pandas as pd

def choose_k(data, max_clus):
    "Helper function to choose the number of clusters for KMeans."
    
    k_tests = np.arange(2, max_clus)
    wcss = []
    silhouettes_ = []
    for k in tqdm.tqdm(k_tests): 
        clustering = KMeans(k, random_state= 42).fit(data)
        wcss.append(clustering.inertia_)
        silhouettes_.append(silhouette_score(data, clustering.predict(data)))
    
    
    return wcss, silhouettes_


def fisher_enrichment_test(df_annot, annotation, cluster, clus_col_name = 'cluster_labels'): 
    """
    Returns a repor dataframe with the top 5 enriched functions
    for a given subset of data. This function is especially suited
    for statistical enrichment tests after clustering. 
    
    Params 
    ------
    df_annot (pd.DataFrame)
        Annotated dataframe containing the 'annotation' column 
        and a 'clus_col_name' column. 
    
    annotation (str)
        Annotation to make the enrichment test on. In the case 
        of gene set enrichment this could be a Gene Ontology 
        or COG annotation. 
    
    cluster (int or str)
        Cluster (or in general group of data points) to test.
    
    col_clus_name (str)
        Name of the cluster column in the df_annot dataframe.
    
    Returns 
    -------
    df_report (pd.DataFrame)
        Report dataframe with pvalues and annotation names. 
    
    """
    # Get subset of completely annotated genes 
    df_test = df_annot[pd.notnull(df_annot[annotation])]

    # Number of genes with valid annotation 
    M = df_test.shape[0]

    # Extract data for given cluster
    df_clus = df_test[df_test[clus_col_name].values == cluster]

    # Get top 5 categories to test
    cats = df_clus[annotation].value_counts().head().index.tolist()
    
    # Number of genes in the cluster (sample size)
    N = df_clus.shape[0]
    
    # Initialize pvalue array 
    pvals = np.empty(len(cats))
    
    # Loop through the top categories
    for i, cat in enumerate(cats): 
        
        df_cat = df_test[df_test[annotation].values == cat]
        
        # Total number of genes that map to given category (total number of white balls)
        n = df_cat.shape[0]
        
        # Number of genes inside cluster that map to given category (number of white balls in sample)
        x = df_clus[df_clus[annotation].values == cat].shape[0]
        
        # Sweep through the probabilities from x to n 
        pmfs = st.hypergeom.pmf(k = np.arange(x, n + 1), N = N, n = n, M = M)
        
        # Compute pval
        pvals[i] = pmfs.sum()
    
    # Save results
    df_report = pd.DataFrame(
        {'categories': cats, 'pval': pvals}
    )

    df_report['cluster'] = cluster
    df_report['annot'] = annotation

    return df_report 

def shuffle_genes(df_not_shuffled):
    ### shuffling of genes with respect to UMAP coordinates:
    cols_genes = df_not_shuffled.columns[:5]
    cols_umap = df_not_shuffled.columns[5:]

    df_genes = df_not_shuffled[cols_genes].copy()
    df_umap = df_not_shuffled[cols_umap].copy()

    df_genes_shuffle = df_genes.reindex(np.random.permutation(df_genes.index)).reset_index(drop = True)
    df_shuffle = pd.concat([df_genes_shuffle, df_umap], axis=1)
    
    return df_shuffle