# mtb_tn_db
*M. tuberculosis* transposon (Tn) screen database 

code subfolders: 

* **data_wrangling:**
  * *data_wrangling.ipynb*: reads in and organizes data from *Tn_datasets* folder into a binary essentiality matrix. 
  * *data_wrangling_qvals_log2FC.ipynb*: reads in and organizes data from *Tn_datasets* folder into a q-value and a log2FC matrix.
  * *Tn_data_wrangling.py*: several functions used by the notebooks above are here. 
  
* **functional_prediction:**
  * *distance_null_hypothesis.ipynb*: exploratory notebook. Computes distributions of Hamming distances between the essentiality profiles / vectors of different subsets of genes. 
  
* **visualizations:**
  * *Tn_library_stats.ipynb*: explores basic statistics of the Tn-matrix. 
  * *counting_unknown_essentials.ipynb*: exploratory notebook. Makes different visualizations of unknown/essential genes.  
  * *visualizing_similarity.ipynb*: exploratory notebook. Makes tSNE-plots using binary essentiality matrix. 

data: 

* **Tn_datasets**: 
  * all the supplementary files I got from publications are stored in this folder.  
  * also contains the genetic interaction datasets from FLUTE knock-out strains. 
* files_and_columns.csv: 
  * specifies for each file, the column name
* files_and_columns_set2.csv
  * specifies for each file, the column name as well as if the file has q-value and log2FC data. 
* column_order.xlsx
  * specifies the order of the columns in the binary matrix.  
* Tn_library_DB.xlsx
  * the current version of the binary essentiality matrix. 
* Tn_library_DB_qval_log2FC.xlsx
  * the current version of the q-value + log2FC matrix. 



