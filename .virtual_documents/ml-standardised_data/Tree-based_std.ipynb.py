import pandas as pd
import pathlib
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.metrics import confusion_matrix
from sklearn.metrics import balanced_accuracy_score, accuracy_score
from sklearn.preprocessing import normalize
#import xgboost as xgb
from collections import defaultdict
from imblearn.over_sampling import SMOTE
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


lfc_mb_filt=pd.read_csv('../data/standardized_data/cleaned_ML/lfc_mb_filt.csv')
lfc_mb_filt.head()


value_cols = [col for col in lfc_mb_filt.columns if col not in ['Rv_ID', 'Functional_Category']]


def rf_vanilla(mat, cols, norm_method, name):
    X=mat[cols].values
    y=mat['Functional_Category'].values
    n_classes=mat['Functional_Category'].nunique()
    #print(X,y)
    accuracy=[]
    confusion=np.zeros((n_classes, n_classes))
    skf = StratifiedKFold(n_splits=3, shuffle=True)
    for train_index, test_index in skf.split(X, y):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        if norm_method=='SMOTE':
            X_train, y_train=SMOTE().fit_resample(X_train, y_train)
            clf=RandomForestClassifier(n_estimators=100)
        if norm_method=='bal_weights':
            clf=RandomForestClassifier(n_estimators=100, class_weight='balanced')
        clf.fit(X_train,y_train)
        y_hat=clf.predict(X_test)
        confusion+=confusion_matrix(y_test, y_hat)
        accuracy.append(accuracy_score(y_test, y_hat))
    print (accuracy)
    
    confusion=confusion/3
    confusion=normalize(confusion, axis=1, norm='l1')
    col_names=['PE/PPE', 'cell wall and\ncell processes', 'information pathways', 'insertion seqs\nand phages', 'intermediary metabolism\nand respiration', 'lipid metabolism', 'regulatory proteins', 'virulence, detoxification,\nadaptation']
    confusion=confusion=pd.DataFrame(confusion, columns=col_names, index=col_names)
    #confusion=pd.DataFrame(confusion, columns=clf.classes_, index=clf.classes_)
    plt.figure(figsize=(7,7))
    rc={'xtick.labelsize': 14, 'ytick.labelsize': 14, 'axes.labelsize': 14}
    sns.set(rc=rc)
    heat=sns.heatmap(confusion, annot=True, linewidths=.1, fmt='1.2f', square=True)
    heat.set(xlabel='PREDICTED CLASS', ylabel='TRUE CLASS', title=name)
    #return accuracy


rf_vanilla(lfc_mb_filt, value_cols, 'SMOTE', 'Accuracy_log2FC_RF_SMOTE')


rf_vanilla(lfc_mb_filt, value_cols, 'bal_weights', 'Accuracy_log2FC_RF_bal_weights')


lfc_bin_mb_filt=pd.read_csv('../data/standardized_data/cleaned_ML/lfc_bin_mb_filt.csv')
lfc_bin_mb_filt.head()


value_cols_lfc_bin = [col for col in lfc_bin_mb_filt.columns if col not in ['Rv_ID', 'Functional_Category']]


rf_vanilla(lfc_bin_mb_filt, value_cols_lfc_bin, 'SMOTE', 'Accuracy_log2FC_bin_RF_SMOTE')


rf_vanilla(lfc_bin_mb_filt, value_cols_lfc_bin, 'bal_weights', 'Accuracy_log2FC_bin_RF_SMOTE')





umap_10 = pd.read_csv('../data/tests/df_umap_10.csv')
umap = pd.read_csv('../data/tests/df_umap.csv')


umap = umap.drop(columns=['gene_name', 'func_tuberculist', 'COG'])
umap_10 = umap_10.drop(columns=['gene_name', 'func_tuberculist', 'COG'])


lfc_umap = lfc_mb_filt.merge(umap, how='left', on='Rv_ID')
lfc_umap_10 = lfc_mb_filt.merge(umap_10, how='left', on='Rv_ID')


#fill empty umaps with zero
lfc_umap = lfc_umap.fillna(0)
lfc_umap_10 = lfc_umap.fillna(0)


lfc_umap_value_cols = [col for col in lfc_umap.columns if col not in ['Rv_ID', 'Functional_Category']]
rf_vanilla(lfc_umap, lfc_umap_value_cols, 'SMOTE', 'Accuracy_log2FC_umap_RF_SMOTE')


lfc_umap_10_value_cols = [col for col in lfc_umap_10.columns if col not in ['Rv_ID', 'Functional_Category']]
rf_vanilla(lfc_umap_10, lfc_umap_10_value_cols, 'SMOTE', 'Accuracy_log2FC_umap_10_RF_SMOTE')
