# Intrinsically linked lineage-specificity of transposable elements and lncRNAs reshapes transcriptional regulation species- and tissue-specifically

This directory contains the eGRAM program and data related to our work: "Intrinsically linked lineage-specificity of transposable elements and lncRNAs reshapes transcriptional regulation species- and tissue-specifically".

Lineage-specificity of transcriptional regulation by lncRNAs critically determines whether mouse models reliably mimic human diseases. To address this question, we identified human/mouse-specific lncRNAs from GENCODE-annotated human and mouse lncRNAs, predicted their DNA binding domains (DBDs) and binding sites (DBSs), analysed transposable elements (TEs) in DBDs and DBSs, and analysed functional enrichment of target genes. 84%/98% of human/mouse-specific lncRNAs, 61%/95% of their DBDs, and 46%/73% of their DBSs contain TEs almost exclusively originated from simians/rodents, indicating intrinsically linked lineage-specificity of TEs, lncRNAs, and lncRNAs’ DBSs. We then revealed how transcriptional regulation is lineage-specifically rewired by co-opted lncRNAs and DBSs by analysing distributions of target genes in signalling pathways and expression of target genes in multiple tissues in humans and mice. Transcriptional regulation is greatly rewired species-specifically and tissue-specifically. We further analysed transcriptomic data of Alzheimer's disease and tumours from human patients and mouse models, with results supporting the above conclusions. Our results reveal the intrinsically linked lineage-specificity of transposable elements, lncRNAs, and transcriptional regulation, provide data and tool for analysing and differentiating transcriptional regulation in humans and mice, and suggest that many evolutionary novelties may be destined to be lineage-specific.

The eGRAM program identifies gene modules comprising co-expressed genes and their regulatory lncRNAs based on lncRNA/DNA bindings and gene expression correlations.

# Requirements
1. **Python**: >=3.7.0

2. **numpy**: >=1.21.6

3. **pandas**: >=1.3.5

4. **scipy**: >=1.7.3

5. **OS**: the eGRAM code has been tested on Linux system.

# Data
1. **data**  --  the TPM-normalized gene expression matrixes derived from patients with Alzheimer's disease (in the lateral temporal lobe, N=12) and mouse models of Alzheimer's disease (in the hippocampus, N=14).

2. **HS_lncRNA_DBS_matrix**  --  the DNA binding matrix of human-specific lncRNAs.

3. **MS_lncRNA_DBS_matrix**  --  the DNA binding matrix of mouse-specific lncRNAs.

4. **PathwayAnnotation**  --  the KEGG pathway annotation of human and mouse.

# Usage
Here is a command line to run the eGRAM program:

```
python eGRAM.py --t1 100 --t2 60 --m 5 --c 0.5 --f1 bindingMatrix --f2 expressionMatrix --f3 KeggCategory --s human --o output
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

# Bug reports
Please send comments and bug reports to JL.linjie@outlook.com.

# Related website
To obtain details about lncRNA/DNA binding prediction using LongTarget, please go to our website http://www.gaemons.net/LongTarget.

# Citation
Lin et al. Intrinsically linked lineage-specificity of transposable elements and lncRNAs reshapes transcriptional regulation species- and tissue-specifically.
