# Unravelling intrinsically linked lineage-specificity of transposable elements, lncRNA genes, and transcriptional regulation

This directory contains the eGRAM program and data related to our work: "Unravelling intrinsically linked lineage-specificity of transposable elements, lncRNA genes, and transcriptional regulation".

The eGRAM program is designed to identify gene modules comprising co-expressed genes and their regulatory lncRNAs based on both lncRNA/DNA binding interactions and gene expression correlations. **The program is available in two versions: eGRAM.py and eGRAMv3R1.py**. The eGRAM.py version is tailored for the analysis of bulk RNA-seq data, while the eGRAMv3R1.py version is specifically developed to handle scRNA-seq data.

# <ins>**eGRAM.py for analyzing bulk RNA-seq data**</ins>
The eGRAM.py (v1) program is a specialized tool designed to identify gene modules that consist of co-expressed genes and their associated regulatory lncRNAs. By integrating both lncRNA/DNA binding interactions and gene expression correlations, eGRAM.py offers a comprehensive approach to analyzing transcriptional regulation. Key features of eGRAM.py include its ability to process bulk RNA-seq data, making it suitable for studies involving tissue-specific or disease-related transcriptional regulation. This version performs an exhaustive search of all possible combinations of transcriptional regulatiors, resulting in increased computational time with large datasets. 

# Requirements
1. **Python**: >=3.7.0

2. **numpy**: >=1.21.6

3. **pandas**: >=1.3.5

4. **scipy**: >=1.7.3

5. **OS**: the eGRAM.py code has been tested on Linux system.

# Data
1. **human_AD_TPM**, **mouse_AD_TPM**  --  the TPM-normalized gene expression matrixes derived from patients with Alzheimer's disease (in the lateral temporal lobe, N=12) and mouse models of Alzheimer's disease (in the hippocampus, N=14).

2. **HS_lncRNA_DBS_matrix**  --  the DNA binding matrix of human-specific lncRNAs.

3. **MS_lncRNA_DBS_matrix**  --  the DNA binding matrix of mouse-specific lncRNAs.

4. **PathwayAnnotation**  --  the KEGG and wiki pathway annotation of human and mouse.

# Usage
Here is a command line to run the eGRAM.py program:

```
python eGRAM.py --t1 100 --t2 60 --m 5 --c 0.5 --f1 HS_lncRNA_DBS_matrix --f2 data/human_AD_TPM --f3 KeggCategory --s human --o output
```

# Help information
Here is a brief explanation of the command line arguments:

```
Options   Parameters      Functions
t1   Binding affinity threshold  An integer, indicate the threshold to determine a strong binding.
t2   Binding affinity threshold  An integer, indicate the threshold to determine a qualified binding.
m    Module size                 An integer, indicate the minimum size of modules.
c    Correlation threshold       A floating point, indicate the minimum correlation coefficient to identify correlated genes.
f1   Binding affinity matrix     A string, indicate the file of binding affinity matrix.
f2   Gene expression matrix      A string, indicate the file of gene expression matrix.
f3   Kegg pathway category       A string, indicate the file of Kegg pathway category.
s    Species                     A string, indicate the species to analyze (e.g., human or mouse).
o    Output                      A string, indicate the output file name.
```

# <ins>**eGRAMv2**</ins>
eGRAMv2 is the basic version and handles RNA-seq data. Employing heuristic search, v2 achieves significantly faster computational speeds than v1 on large datasets. The inputs include (a) an RNA-seq dataset, (b) a lncRNA-DBS file containing lncRNAs and predicted targets, (c) a TF-DBS file containing TFs and predicted targets. If the user just examines transcriptional regulation by lncRNAs, he/she can slightly revise the code to make (c) as optional. In addition, the KEGG and WikiPathways databases that hold human and mouse pathways are the built-in inputs.

Main parameters include (a) a lncRNA-DBS threshold, (b) a TF-DBS threshold, (c) a module size threshold, (d) species (1=human, 2=mouse), which determines using human or mouse KEGG/WikiPathways pathways, (e) a correlation threshold for lncRNA/target and lncRNA/lncRNA, (f) a correlation threshold for TF/target and TF/TF, (g) Pearson/Spearman that determines computing Pearson coefficient or Spearman coefficient as the measure of correlation, (h) a FDR threshold for determining significance of pathway enrichment (default 0.01), (i) a fold name for outputting files of results. All default parameter values are given in the code.

eGRAMv2 performs the following steps.
 
step 1 - read in user files and command line arguments

step 2 - read in kegg files and wikipathway files

step 3 - computing lncRNA-gene correlations, lncRNA-lncRNA correlations, and upon which, lncRNAsets

step 4 - computing TF-gene correlations, TF-TF correlations, and upon which, TFsets

step 5 - generate shared lncRNA target sets upon df_lncRNA_gene_DBS and df_lncRNA_gene_corr

step 6 - generate shared TF target sets upon df_TF_gene_DBS and df_TF_gene_corr

step 7 - identify typically co-regulated modules and independently-regulated modules

step 8 - perform pathway enrichment analysis for module genes in df_path_module

step 9 - write results to files

step 10 - generate files for cytoscape

eGRAMv2 can be run in a terminal or within an IDE (integrated development environment) such as Pycharm.  

# Requirements
1. **Python**: >=3.7.0

2. **numpy**: >=1.21.6

3. **pandas**: >=1.3.5

4. **scipy**: >=1.7.3

5. **OS**: the eGRAMv2R1 code has been tested on Linux system.

# Data
1. **2024May-DEG-exp-A549-2513.csv**  --  the gene expression matrix.

2. **2024May-lncRNA-DBS-A549-2513.csv**  --  the DNA binding matrix (i.e. target genes) of lncRNAs.

3. **2024May-TF-DBS-A549-2513.csv**  --  the DNA binding matrix (i.e. target genes) of TFs.

4. **KEGG**  --  the KEGG pathway annotation of human and mouse.

5. **WikiPathways**  --  the WikiPathways annotation of human and mouse. The compressed package needs to be decompressed before running the program.

# Usage
Here is a command line to run the eGRAMv2R1 program:

```
'Example: python eGRAMv2R1.py --exp 2024May-DEG-exp-A549-2513.csv --lncDBS 2024May-lncRNA-DBS-A549-2513.csv --tfDBS 2024May-TF-DBS-A549-2513.csv --species 1  --lncCutoff 100 --tfCutoff 8 --lncCorr 0.5 --tfCorr 0.5 --moduleSize 50 --corr Pearson --remove_redundancy y --fdr 0.05  --out myout'
```

# Help information
Here is a brief explanation of the command line arguments:

```
Parameters      Functions
--exp gene_exp_file is required; if genetype does not have TF, all TF-related functions are not performed. 
--lncDBS     lncRNA_DBS_file is required; integer or 0. 
--tfDBS      TF_DBS_file is not a must; required if gene_exp_file contains TF.
--species    1=human, 2=mouse.
--lncCutoff  the DBS threshold for determining lncRNAs' targets (defult=100).
--tfCutoff   the DBS threshold for determining TFs's targets (default=10).
--lncCorr    the correlation threshold for lncRNA/target and lncRNA/lncRNA (default=0.5).
--tfCorr     the correlation threshold for TF/target and TF/TF (default=0.5).
--moduleSize the minimum gene number of a module (default=50).
--corr       Pearson/Spearman (default=pearson).
--remove_redundancy    y/n, whether or not merge modules with the same gene sets (default=y)
--fdr        the FDR value to determine significance of pathway enrichment (default 0.01).
--out        output file name.
```


# <ins>**eGRAMv3 for analyzing scRNA-seq data**</ins>
eGRAMv3 is specifically designed for handling scRNA-seq data. The minimal inputs include (a) an scRNA-seq dataset, (b) a lncRNA-DBS file containing lncRNAs and predicted targets. A TF-DBS file containing TFs and predicted targets is an optional. As The KEGG and WikiPathways databases are built-in components. 

Parameters include (a) a lncRNA-DBS threshold, (b) a TF-DBS threshold (optional), (c) a module size threshold, (d) a zero proportion threshold (for filtering genes with too high dropouts), (e) a FDR threshold for determining significance of pathway enrichment (default 0.01), (f) a fold name for outputting files of results, (g) y/n indicating whether merging sets or not, (h) a threshold for merging sets (when the percentage of shared genes between two models is lower than the threshold, stopping merging modules), (i) species (1=human, 2=mouse), which determines using human or mouse KEGG/WikiPathways pathways, (j) the alpha parameter, which is important for controlling the MIC algorithm. 

eGRAMv3 performs the following steps. (a) For each lncRNA, identity its partners upon expression correlation and generate its lncRNA set. (b) For each lncRNA set, compute correlations between lncRNAs therein and all genes and generate a target set of this lncRNA set. (c) Combine the targeting and correlation relationships and generate target sets for all lncRNA sets. (d) Merge lncRNA sets that share the same target set into a regulator set and build a module consisting of the regulatory set and the target set. (e) In a module, a lncRNA is a activator if it has positive correlations with >80% of the targets, or a repressor if it has negative correlations with >80% of the targets. (f) Determine whether a module is significantly enriched for genes in pathways in the KEGG and Wikipathway databases. (g) Modules can be further merged if their target sets share substantial components. In each round of merging, two modules sharing the highest percentage of targets are merged; and if >=2 modules have the same percentage, the two modules sharing the highest percentage of regulators are merged. Shared components in the regulatory set and/or target set across modules connect modules into a network.

scRNA-seq data are notorious for abundant dropouts, and various algorithms and tools have been developed either for pro-processing dropouts or for downstream analysis. We note that dropouts do not impair the ability of eGRAM to identify regulatory modules, but they greatly influence the computation of correlation. Lacking control samples is another feature that makes scRNA-seq analysis differ from RNA-seq analysis. eGRAMv3 can process normalized counts of scRNA-seq data with dropouts. It uses the Maximal In-formation Coefficient (MIC) algorithm to compute correlations between genes. MIC belongs to the maximal information-based nonparametric exploration (MINE) class of statistics, can capture a wide range of correlations by discretizing the data and finding the grid that maximizes the mutual information between variables, has relatively equal power across different types of relationships, and can tolerate dropouts if sample size is large. 

MIC returns two estimators - mic (the maximal information coefficient) and tic (total information coefficient). Since a large mic or tic value alone does not necessarily indicate correlation because they are influenced by sample size and grid resolution, the original algorithm uses permutation test to examine the significance of mic and tic values. For each pair of variables, the test generates a null distribution by shuffling the data, recalculates mic and tic many times, and compares the observed mic against this null distri-bution. This method does not fit scRNA-seq data because of too many variables. Since the majority of gene pairs are uncorrelated, the bulk of the mic values can be used to estimate the null distribu-tion, and mic values come from a mixture of null (background) and non-null (correlated) distribu-tions. Our method explores adaptive thresholding via null distribution modeling, which extracts upper triangle values from the mic matrix, fits a Gaussian mixture model with two components, and get the threshold as “mean + 3 standard deviation” of background component, and identifies signif-icant pairs upon the threshold. There are multiple ways to use the mic and tic values: using mic alone, using mic AND tic, and using mic OR tic. The user can make choice upon testing running.

Since many scRNA-seq datasets contain many genes and cells, a parallel version of eGRAMv3 is developed to explore all cores in a CPU. When running this version, the compute may respond very slowly to keyboard and mouse input. If the user does not want to use the version, he/she can slightly revise the code, chooing pstats() instead of par_pstats(). 

The alpha parameter controls the for-loop that finds the grid that maximizes the mutual information between variables, it thus controls time consumption of the MIC algorithm. In addition the parallel version, adaptive alpha values can be set depending on gene and cell numbers in a dataset. 

# Requirements
1. **Python**: >=3.7.0

2. **numpy**: >=1.21.6

3. **pandas**: >=1.3.5

4. **scipy**: >=1.7.3

5. **sklearn**: >=1.0.2

6. **OS**: the eGRAM code has been tested on Linux and Windows system.

# Data
1. **M_SPG_2025.csv**  --  the scRNA-seq data of mouse SPG cells.

2. **Hmarker_HlncRNA_DBS.csv**  --  the DNA binding matrix of human lncRNAs.

3. **Mmarker_MlncRNA_DBS.csv**  --  the DNA binding matrix of mouse lncRNAs.

4. **PathwayAnnotation**  --  the KEGG and wiki pathway annotation of human and mouse.

# Usage
Here is a command line to run the eGRAMv3 program:

```
python eGRAMv3R1.py --exp data/M_SPG_2025.csv --lncDBS data/Mmarker_MlncRNA_DBS.csv --lncCutoff 36 --moduleSize 50 --fdr 0.01 --zeroP 0.85 --merge y --moduleOverlap 0.98 --outFold M-SG-1 --species 2 --alpha 0.0 --outFold out
```

# Help information
Here is a brief explanation of the command line arguments:

```
Parameters:
-exp             scRNA-seq file. 
-lncDBS          lncRNA DBS file. 
-tfDBS           TF DBS file.
-lncCutoff       lncRNA DBS binding affinity threshold (lncRNA DBSs < this threshold will be ignored).
-tfCutoff        TF DBS binding affinity threshold (TF DBSs < this threshold will be ignored).
-moduleSize      module size threshold (modules < this threshold will be ignored).
-fdr             FDR for pathway enrichment analysis (pathways whose significance > this threshold will be ignored).
-zeroP           zero proportion threshold (genes whose zero proportion > this threshold will be ignored).
-merge           whether merge modules with highly redundant genes (y=merge, n=not). 
-moduleOverlap   module overlapping threshold (stop merging when overlap between modules < this threshold).
-outFold         output fold.
-species         1 = human; 2 = mouse
-alpha           XXX = fixed alpha; 0.0 = an alpha depending on the detected sample size.  
```

# Bug reports
Please send comments and bug reports to JL.linjie@outlook.com or zhuhao@smu.edu.cn.

# Related website
To obtain details about lncRNA/DNA binding prediction using LongTarget, please go to our website http://www.gaemons.net/LongTarget.

# Citation
Lin et al. Unravelling intrinsically linked lineage-specificity of transposable elements, lncRNA genes, and transcriptional regulation.
