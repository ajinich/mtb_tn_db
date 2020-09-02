def lasso_confusion(mat, cols, C, norm_method, name, figsize=(7, 7)):
    X = mat[cols].values
    y = mat['Functional_Category'].values
    n_classes = mat['Functional_Category'].nunique()
    accuracy = []
    n_splits = 3
    confusion = np.zeros((n_classes, n_classes))
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    for train_index, test_index in skf.split(X, y):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        if norm_method == 'SMOTE':
            X_train, y_train = SMOTE().fit_resample(X_train, y_train)
            clf = LogisticRegression(
                penalty='l1', solver='liblinear', multi_class='ovr', C=C, random_state=42)
        elif norm_method == 'bal_weights':
            clf = LogisticRegression(
                penalty='l1', solver='liblinear', multi_class='ovr', C=C, class_weight='balanced', random_state=42)
        elif norm_method == None:
            clf = LogisticRegression(
                penalty='l1', solver='liblinear', multi_class='ovr', C=C, random_state=42)
        clf.fit(X_train, y_train)
        y_hat = clf.predict(X_test)
        confusion += confusion_matrix(y_test, y_hat)
        accuracy.append(accuracy_score(y_test, y_hat))
    confusion = normalize(confusion, axis=1, norm='l1')
    print(accuracy)

    # Create confusion matrix
    col_names = ['PE/PPE', 'cell wall and\ncell processes', 'information pathways', 'insertion seqs\nand phages',
                 'intermediary metabolism\nand respiration', 'lipid metabolism', 'regulatory proteins', 'virulence, detoxification,\nadaptation']
    confusion = confusion = pd.DataFrame(
        confusion, columns=col_names, index=col_names)
    #confusion=pd.DataFrame(confusion, columns=clf.classes_, index=clf.classes_)
    plt.figure(figsize=figsize)
    rc = {'xtick.labelsize': 14, 'ytick.labelsize': 14, 'axes.labelsize': 14}
    sns.set(rc=rc)
    heat = sns.heatmap(confusion, annot=True, linewidths=.1,
                       fmt='1.2f', square=True)
    heat.set(xlabel='PREDICTED CLASS', ylabel='TRUE CLASS', title=name)
    return None


def lasso_coefs(mat, cols, C, norm_method):
    X = mat[cols].values
    y = mat['Functional_Category'].values
    # print(X,y)
    if norm_method == 'SMOTE':
        X, y = SMOTE().fit_resample(X, y)
        clf = LogisticRegression(
            penalty='l1', solver='liblinear', multi_class='ovr', C=C, random_state=42)
    if norm_method == 'bal_weights':
        clf = LogisticRegression(
            penalty='l1', solver='liblinear', multi_class='ovr', C=C, class_weight='balanced', random_state=42)
    clf.fit(X, y)
    coefs = clf.coef_
    coefs = np.transpose(coefs)
    coefs = pd.DataFrame(coefs, columns=clf.classes_, index=cols)
    # print(clf.coef_)
    return coefs


def lasso_vanilla(mat, cols, Cs):
    X = mat[cols].values
    y = mat['Functional_Category'].values
    # print(X,y)
    accuracy = {}
    for C in Cs:
        skf = StratifiedKFold(n_splits=3, shuffle=True)
        accuracy_per_fold = []
        for train_index, test_index in skf.split(X, y):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]
            X_train, y_train = SMOTE().fit_resample(X_train, y_train)
            clf = LogisticRegression(
                penalty='l1', solver='liblinear', multi_class='ovr', C=C)
            clf.fit(X_train, y_train)
            y_hat = clf.predict(X_test)
            accuracy_per_fold.append(balanced_accuracy_score(y_test, y_hat))
        accuracy[C] = sum(accuracy_per_fold)/len(accuracy_per_fold)
    print(accuracy)
    # return accuracy


get_ipython().run_line_magic("reset", " -f")


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from imblearn.over_sampling import SMOTE
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.model_selection import StratifiedKFold

get_ipython().run_line_magic("matplotlib", " inline")


lfc_mb_filt = pd.read_csv(
    '../data/standardized_data/cleaned_ML/lfc_mb_filt.csv')
qval_mb_filt = pd.read_csv(
    '../data/standardized_data/cleaned_ML/qval_mb_filt.csv')


# replaced values of qval lower than 0.05 with 0.05. Then divided log2FC by qval^beta
def lasso_inv_p(betas, Cs, lfc, qval, value_cols):
    X = lfc[value_cols].values
    y = lfc['Functional_Category'].values
    # replace values lower than 0.05 with 0.05
    qval = qval[value_cols].apply(lambda x: np.where(x < 0.05, 0.05, x))
    qval = qval[value_cols].values
    # print(X,y)
    accuracy = pd.DataFrame(index=betas, columns=Cs, dtype='float')
    for b in betas:
        X = X/np.power(qval, b)
        for C in Cs:
            # print(b)
            skf = StratifiedKFold(n_splits=3, shuffle=True)
            accuracy_per_fold = []
            for train_index, test_index in skf.split(X, y):
                X_train, X_test = X[train_index], X[test_index]
                y_train, y_test = y[train_index], y[test_index]
                X_train, y_train = SMOTE().fit_resample(X_train, y_train)
                clf = LogisticRegression(
                    penalty='l1', solver='liblinear', multi_class='ovr', C=C)
                clf.fit(X_train, y_train)
                y_hat = clf.predict(X_test)
                accuracy_per_fold.append(accuracy_score(y_test, y_hat))
            accuracy.loc[b, C] = sum(accuracy_per_fold)/len(accuracy_per_fold)
    plt.figure()
    heat = sns.heatmap(accuracy, annot=True, linewidths=.1, fmt='1.2f')
    heat.set(xlabel='Cs', ylabel='betas',
             title='Accuracy with LR on log2FC w pval')
    return None


betas = [0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4]
Cs = [0.005, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100]





(lfc_mb_filt.columns == qval_mb_filt.columns).all()


lfc_mb_filt.head()


value_cols = [col for col in lfc_mb_filt.columns if col not in [
    'Rv_ID', 'Functional_Category']]


lasso_inv_p(betas, Cs, lfc_mb_filt, qval_mb_filt, value_cols)


def lasso_min_p(lfc, qval, C, value_cols):
    X = lfc[value_cols].values
    y = lfc['Functional_Category'].values
    # replace values lower than 0.05 with 0.05
    qval = qval[value_cols].apply(lambda x: np.where(x < 0.05, 0.05, x))
    qval = qval[value_cols].values
    # print(X,y)
    accuracy = []
    X = X*(1-qval)
    skf = StratifiedKFold(n_splits=3, shuffle=True)
    for train_index, test_index in skf.split(X, y):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        X_train, y_train = SMOTE().fit_resample(X_train, y_train)
        clf = LogisticRegression(
            penalty='l1', solver='liblinear', multi_class='ovr', C=C)
        clf.fit(X_train, y_train)
        y_hat = clf.predict(X_test)
        accuracy.append(accuracy_score(y_test, y_hat))
    print(accuracy)
    # return accuracy


lasso_min_p(lfc_mb_filt, qval_mb_filt, 10, value_cols)





get_ipython().run_line_magic("reset", " -f")


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from imblearn.over_sampling import SMOTE
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import normalize

get_ipython().run_line_magic("matplotlib", " inline")


lfc_mb_filt = pd.read_csv(
    '../data/standardized_data/cleaned_ML/lfc_mb_filt.csv')


value_cols = [col for col in lfc_mb_filt.columns if col not in [
    'Rv_ID', 'Functional_Category']]


lasso_confusion(lfc_mb_filt, value_cols, C=10,
                norm_method='SMOTE', name='Accuracy_log2FC_LR_SMOTE_C=10')


lasso_confusion(lfc_mb_filt, value_cols, C=10, norm_method='bal_weights',
                name='Accuracy_log2FC_LR_balweights_C=10')


lasso_confusion(lfc_mb_filt, value_cols, C=10, norm_method=None,
                name='Accuracy_log2FC_LR_NoNorm_C=10')


lfc_mb_filt.Functional_Category.value_counts()


lasso_coefs(lfc_mb_filt, value_cols, 10, 'SMOTE').to_csv(
    'results/coefs_mb_lfc_SMOTE_C10.csv')





get_ipython().run_line_magic("reset", " -f")


import pandas as pd
import pathlib
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import confusion_matrix, accuracy_score, balanced_accuracy_score
from sklearn.preprocessing import normalize
from collections import defaultdict
from imblearn.over_sampling import SMOTE
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic("matplotlib", " inline")


bin_mb_filt = pd.read_csv(
    '../data/standardized_data/cleaned_ML/bin_mb_filt.csv')


binary_data_cols = [col for col in bin_mb_filt.columns if col not in [
    'Rv_ID', 'Functional_Category']]


lasso_vanilla(bin_mb_filt, binary_data_cols, Cs=[0.005, 0.05, 0.1, 1, 10, 50])


lasso_confusion(bin_mb_filt, binary_data_cols, C=10,
                norm_method='SMOTE', name='Accuracy_bin_LR_SMOTE_C=10')


lasso_confusion(bin_mb_filt, binary_data_cols, C=10,
                norm_method='bal_weights', name='Accuracy_bin_LR_balweights_C=10')





get_ipython().run_line_magic("reset", " -f")


import pandas as pd
import pathlib
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.metrics import confusion_matrix, balanced_accuracy_score, accuracy_score
from sklearn.preprocessing import normalize
from collections import defaultdict
from imblearn.over_sampling import SMOTE
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic("matplotlib", " inline")


lfc_mb_filt = pd.read_csv(
    '../data/standardized_data/cleaned_ML/lfc_mb_filt.csv')
bin_mb_filt = pd.read_csv(
    '../data/standardized_data/cleaned_ML/bin_mb_filt.csv')
assert (lfc_mb_filt.columns == bin_mb_filt.columns).all()


value_cols = [col for col in lfc_mb_filt.columns if col not in [
    'Rv_ID', 'Functional_Category']]


bin_mb_filt = bin_mb_filt.rename(
    columns={col: col + '_bin' for col in value_cols})


lfc_bin_mb = pd.merge(lfc_mb_filt, bin_mb_filt, how='left', on='Rv_ID').rename(columns={
    'Functional_Category_y': 'Functional_Category'}).drop(columns=['Functional_Category_x'])


lfc_bin_mb.to_csv(
    '../data/standardized_data/cleaned_ML/lfc_bin_mb_filt.csv', index=False)


lfc_bin_cols = [col for col in lfc_bin_mb.columns if col not in [
    'Rv_ID', 'Functional_Category']]


lasso_vanilla(lfc_bin_mb, lfc_bin_cols, Cs=[0.005, 0.05, 0.1, 1, 10, 50])


lasso_confusion(lfc_bin_mb, lfc_bin_cols, C=10,
                norm_method='SMOTE', name='Accuracy_lfc_bin_LR_SMOTE_C=10')


get_ipython().run_line_magic("reset", " -f")


import pandas as pd
import pathlib
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.metrics import confusion_matrix, balanced_accuracy_score, accuracy_score
from sklearn.preprocessing import normalize
from collections import defaultdict
from imblearn.over_sampling import SMOTE
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic("matplotlib", " inline")


umap_10 = pd.read_csv('../data/tests/df_umap_10.csv')
umap = pd.read_csv('../data/tests/df_umap.csv')


umap.head()


umap = umap.drop(columns=['gene_name', 'func_tuberculist', 'COG'])
umap_10 = umap_10.drop(columns=['gene_name', 'func_tuberculist', 'COG'])


lfc_mb_filt = pd.read_csv(
    '../data/standardized_data/cleaned_ML/lfc_mb_filt.csv')
lfc_mb_filt.head()


lfc_umap = lfc_mb_filt.merge(umap, how='left', on='Rv_ID')
lfc_umap_10 = lfc_mb_filt.merge(umap_10, how='left', on='Rv_ID')


lfc_umap.isna().sum()


#fill empty umaps with zero
lfc_umap = lfc_umap.fillna(0)
lfc_umap_10 = lfc_umap.fillna(0)



lfc_umap_value_cols = [col for col in lfc_umap.columns if col not in ['Rv_ID', 'Functional_Category']]
lasso_confusion(lfc_umap, lfc_umap_value_cols, C=10, norm_method='SMOTE',
                name='Accuracy_log2FC_umap_LR_SMOTE_C=10')


lfc_umap_10_value_cols = [col for col in lfc_umap_10.columns if col not in ['Rv_ID', 'Functional_Category']]
lasso_confusion(lfc_umap_10, lfc_umap_10_value_cols, C=10, norm_method='SMOTE',
                name='Accuracy_log2FC_umap_10_LR_SMOTE_C=10')



