# Transposable elements rewire transcriptional regulation by mediating co-option between human/mouse-specific lncRNAs and target genes

This directory contains the eGRAM program and data related to our work: "Transposable elements rewire transcriptional regulation by mediating co-option between human/mouse-specific lncRNAs and target genes".

In this work, we identified 66 and 212 human- and mouse-specific lncRNAs from GENCODE-annotated human and mouse lncRNAs, predicted their DNA binding domains (DBDs) and binding sites (DBSs), identified TEs in these DBDs and DBSs, and analyzed the impacts of TEs and human/mouse-specific lncRNAs on transcriptional regulation. Transcriptional analyses indicate that simian/rodent TEs rewire transcriptional regulation human- and mouse-specifically by mediating the co-option between human/mouse-specific lncRNAs and their target genes.

The eGRAM program identifies gene modules comprising co-expressed genes and their regulatory lncRNAs based on lncRNAs/DNA bindings and gene expression correlations.

# Requirements
1. **Python**: >=3.7.0

2. **numpy**: >=1.21.6

3. **pandas**: >=1.3.5

4. **scipy**: >=1.7.3

5. **OS**: the eGRAM code has been tested on Linux system.

# Data
1. **logTPM_matrix**  --  the log-transformed gene expression matrixes of the 13 human GTEx brain regions and the mouse brain.

2. **HS_lncRNA_DBS_matrix**  --  the DNA binding matrix of human-specific lncRNAs.

3. **MS_lncRNA_DBS_matrix**  --  the DNA binding matrix of mouse-specific lncRNAs.

4. **PathwayAnnotation**  --  the KEGG pathway annotation of human and mouse.

# Bug reports
Please send comments and bug reports to JL.linjie@outlook.com.

# Related website
To obtain details about lncRNA/DNA binding prediction using LongTarget, please go to our website http://www.gaemons.net/LongTarget.

# Citation
Lin et al. Transposable elements rewire transcriptional regulation by mediating co-option between human/mouse-specific lncRNAs and target genes.
