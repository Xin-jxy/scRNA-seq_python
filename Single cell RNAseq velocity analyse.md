## Single cell RNAseq velocity analyse

This README is an indication of different script's utilization, which contains a script of pre-processing of raw dataset with scanpy`single-cell-preprocessing.py `; a script of analyse RNAseq for the distinction of different cluster with scanpy`single-cell-analysis.py`, which could also be done with scVelo; a script of velocity analysis`sc-velocity.py`, which add raw spliced, unspliced and ambiguous data to raw dataset and pre-processing the dataset with same parameter like `single-cell-preprocessing.py `,and then the velocity analysis, the main script of this internship,  `parametre_sc_velocity.py`is the script who can manipulate `parser_of_data.py` is an exclusif script for the annotation of new cluster in terms of clusters annotated primairly with D2 and B6 linear and their barcodes associated. 


### **single-cell-preprocessing.py **

exemple of code:

`python3 single-cell-preprocessing.py  path/to/dataset/repertory/BM/Gene/raw/` 

if you want to generate the list of barcode associated or the result of filteration of dataset:

`python3 single-cell-preprocessing.py barcodes.csv/result_BM.h5ad`

if you want them both:

`python3 single-cell-preprocessing.py result_BM.h5ad barcodes.csv `


### **single-cell-analysis.py**

exemple of code:

`python3 single-cell-analysis.py linear_B6.h5ad`

In order to concatenate 3 samples of a linear:

`python3 single-cell-analysis.py sample_BW.h5ad sample_BM.h5ad sample_B6.h5ad linear_B6.h5ad`

In order to concatenate 2 linears or 2 samples:

`python3 single-cell-analysis.py linear_B6.h5ad linear_D2.h5ad linear_B6_D2.h5ad`
or
`python3 single-cell-analysis.py sample_BW.h5ad sample_BM.h5ad linear_BM.h5ad`

### **sc-velocity.py and parametre_sc_velocity.py**

NB: `parametre_sc_velocity.py` may have this import module problem, use `sc_velocity_back_up.py` if it's the case.

exemple of code:

1. Generate filtered dataset interms of raw dataset(skip this step if you've already dataset filtered by Seurat or any other ways)

`python3 parametre_sc_velocity.py path/to/dataset/repertory/BM/Gene/raw/ file/to/spliced.mtx file/to/unspliced.mtx file/to/ambiguous.mtx new_dataset.h5ad`

2. Input different dataset filtered by following different indication

`python3 parametre_sc_velocity.py dataset_filtered.h5ad`

**NB:The parametre in the plot needed be changed in terms of different demands, here the plot is in terms of louvain calculate(color or groupby in the function), I've not tired cuz too big in my computer and still have error in Linux,which we never fixe (sorry to throw this ball on you, william), also the annotation is based on the filtered dataset, so if you start a new analysis, the ancien annotation doesn't work anymore.**

3. There is the possibility for generating barcodes associated for using in `parser_of_data.py`

`python3 parametre_sc_velocity.py dataset_filtered.h5ad barcodes.csv`

NB: Also if we want to rearrange the cluster, there is always a function of cellrank available. 

### **`parser_of_data.py`**

This script allows tell you the new annotation which have already annotated at the first time in scanpy, please verify your annotation before you use this script. All the basic annotation is writen in `single-cell-analysis.py`, the annotation of cluster is writen in `sc-velocity.py` is based on the jugement of this script. 

exemple of code:

`python3 barcodes_of_scanpy_annotated_cluster.csv barcodes_of_scVelo_non_annotated.csv`


