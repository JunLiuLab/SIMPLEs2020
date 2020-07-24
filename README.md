# SIMPLEs2020: Reproducible codes for [SIMPLEs](https://github.com/JunLiuLab/SIMPLEs) experiments

## Data sets used in the [manuscript](https://www.biorxiv.org/content/10.1101/2020.01.13.904649v1?rss=1)
* Can be downloaded from the Zenodo: https://doi.org/10.5281/zenodo.3958371

## Explanations about the mESC data and hECS data in Zenodo.
NOTE: Links below to GitHub might be broken. Please see the corresponding codes and data in the Zenodo.
* mESC (Deng's) data
  * I used this dataset to show that: 
    *  1) SIMPLEs can discover subtypes of cells. TSNE plot of imputed data
    *  2) Imputation result of differential expressed genes. Heatmap

  * Raw data: https://hemberg-lab.github.io/scRNA.seq.datasets/mouse/edev/#deng
  * Preprocessed data:
    * [deng_dat](mESC_deng_dataset/deng_dat.rdat): log-normalized and filtered an outlier cell and genes.  
    * It has 4 objects: celltype_true (6 major cell types), cl2 (10 subtypes), Y2 (preprocessed data), clus(initial clustering results by Kmeans)

  * The script for preprocessing and plotting the results: [plot](mESC_deng_dataset/deng_preprocess_plot.R)
  * The script for different imputation methods: [SCRABBLE_VIPER_SAVER](mESC_deng_dataset/deng_SCRABBLE_VIPER_SAVER.R), [SIMPLE](mESC_deng_dataset/deng_SIMPLE.R) 
  and the [results](mESC_deng_dataset/Deng).

* hECS (Chu's) data
  * Chu’s data has two parts: 1) chu1: different cell types; 2) chu2: time series. Both of them have corresponding bulk RNASeq. Chu’s data has good quality, so I added dropouts to the original data to compare if can recover the truth. 

  * Raw data: downloaded from GEO: GSE75748
  * Preprocessed data: 
    * [Chu_celltype](hESC_chu_dataset/chu_data.rdat): log-normalized and filtered genes for part 1. 
      * [preprocess](hESC_chu_dataset/chu1_celltype_preprocess.R)
      * It has 5 objects: bulk_norm (log-normalized bulk RNASeq), dat_norm (log-normalized scRNASeq), celltype (cell type label for scRNASeq), celltype0 (cell type label for bulk RNASeq), bulk_norm_mean (mean of gene expression for each cell type from bulk RNASeq). The rows (genes) of bulk and scRNASeq data are the same.

    * [Chu_ts](hESC_chu_dataset/chu_data_ts.rdat): log-normalized and filtered genes for part 2. 

  * Chu's cell type data (Chu1) based simulation, i.e adding dropouts to original data set: 
    * Rscripts: [SCRABBLE_VIPER_SAVER](hESC_chu_dataset/chu1_SCRABBLE_VIPER_SAVER.R), [SIMPLE_MAGIC_SCIMPUTE](chu1_SIMPLE_MAGIC_SCIMPUTE.R)
  
  * Chu's cell type data: 
    * all cells (Chu_all):
      * Rscript: [SCRABBLE_VIPER_SAVER](hESC_chu_dataset/chu_all_more_methods.R), [MAGIC_SCIMPUTE](hESC_chu_dataset/chu_all.R), [SIMPLES](hESC_chu_dataset/chu_all_new.R)
      * Results: [Others_result/](Others_result)Chu_all_*.RData, [SIMPLEs](SIMPLES_result/chu_all2_0.8_1_10.rdat)
    * only DE and EC cell types: [plot](hESC_chu_dataset/chu_analysis_clean.R)
  
  * Chu's time series data (Chu2): 
    * Rscripts: [SCRABBLE_VIPER_SAVER](hESC_chu_dataset/chu2_SCRABBLE_VIPER_SAVER.R), [MAGIC_SCIMPUTE](hESC_chu_dataset/chu2_clean.R), [SIMPLES](hESC_chu_dataset/chu_ts_new.R)
    * Results: [Others_result/](Others_result)Chu_ts_*_.RData, scimpute_ts*, [SIMPLE-B](SIMPLES_result/chu_ts_bulk_0228-099.rdat)

## Reference
SIMPLEs: a single-cell RNA sequencing imputation strategy preserving gene modules and cell clusters variation
Zhirui Hu, Songpeng Zu, Jun S. Liu, bioRiv 2020
