# Unravelling intrinsically linked lineage-specificity of transposable elements, lncRNA genes, and transcriptional regulation

This directory contains the eGRAM program and data related to our work: "Unravelling intrinsically linked lineage-specificity of transposable elements, lncRNA genes, and transcriptional regulation".

The eGRAM program is designed to identify gene modules comprising co-expressed genes and their regulatory lncRNAs based on both lncRNA/DNA binding interactions and gene expression correlations. The program is available in two versions: eGRAM.py and eGRAMv3R1.py. The eGRAM.py version is tailored for the analysis of bulk RNA-seq data, while the eGRAMv3R1.py version is specifically developed to handle scRNA-seq data.

# <ins>**eGRAM.py for analyzing bulk RNA-seq data**</ins>
The eGRAM.py program is a specialized tool designed to identify gene modules that consist of co-expressed genes and their associated regulatory lncRNAs. By integrating both lncRNA/DNA binding interactions and gene expression correlations, eGRAM.py offers a comprehensive approach to analyzing transcriptional regulation. Key features of eGRAM.py include its ability to process bulk RNA-seq data, making it suitable for studies involving tissue-specific or disease-related transcriptional regulation.

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


# <ins>**eGRAMv3 for analyzing scRNA-seq data**</ins>
eGRAMv3 is an advanced version of the eGRAM program, specifically engineered to handle the complexities of scRNA-seq data. It builds on the core functionality of identifying gene modules that consist of co-expressed genes and their associated regulatory lncRNAs. However, eGRAMv3 introduces several enhancements to accommodate the unique challenges of scRNA-seq data, such as high dropout rates and the need for more nuanced analysis of gene expression at the single-cell level.

Key features of eGRAMv3 include:

1. **Handling Dropout Events**: eGRAMv3 is designed to effectively manage the high proportion of dropout events commonly observed in scRNA-seq data. It ensures that these technical zeros do not disproportionately affect the identification of gene modules and regulatory relationships.

2. **Improved Correlation Analysis**: The program uses the Maximal Information Coefficient (MIC) algorithm to compute correlations between genes. MIC is particularly well-suited for scRNA-seq data as it can capture a wide range of non-linear relationships and is robust to the high noise levels typical of single-cell data.

3. **Flexible Parameter Adjustment**: eGRAMv3 allows for fine-tuning of parameters to suit the specific characteristics of different scRNA-seq datasets. This flexibility enables researchers to optimize the analysis for their particular study, whether it involves a large number of cells, a high degree of cellular heterogeneity, or other complex data structures.

4. **Comprehensive Regulatory Network Construction**: Like its predecessor, eGRAMv3 constructs regulatory networks by integrating lncRNA/DNA binding interactions with gene expression correlations. However, it does so with enhanced precision and specificity, allowing for the identification of subtle regulatory relationships that might be missed in bulk RNA-seq analyses.

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
