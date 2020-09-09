import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
get_ipython().run_line_magic("matplotlib", " inline")


confusion = pd.DataFrame({'CE':[1157, 1205], 'not_CE':[1292, 107952]}, index=['CE', 'not_CE'])


confusion


plt.figure(figsize=(7,7))
rc={'xtick.labelsize': 14, 'ytick.labelsize': 14, 'axes.labelsize': 14}
sns.set(rc=rc)
heat=sns.heatmap(confusion, annot=True, linewidths=.1, fmt='1.2f', square=True)
heat.set(xlabel='SI_data', ylabel='Std_data')
fig = heat.get_figure()
fig.savefig('standardized_data_confusion.png', dpi=500, bbox_inches = "tight")
