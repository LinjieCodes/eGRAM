'''
The parallel code is occasionally unstable, especially under Pycharm in the step-run mode, and reports errors such as
"Process finished with exit code -1073741819 (0xC0000005)" or "BrokenPipeError". Re-run the task when it occurs.
'''

#!/usr/bin/python
import multiprocessing
import itertools
import numpy as np
import pandas as pd
import re
import os
import csv
import sys, getopt
from itertools import combinations
from math import log, sqrt
from sklearn.mixture import GaussianMixture
from scipy.stats import gaussian_kde
from scipy.stats import hypergeom
from scipy.stats import rankdata

########################################################################################################################
# Auxiliary Functions - MIC algorithm
class MINE:
    def __init__(self, alpha, c=1, est="mic_e"):
        """
        Python implementation of MINE object
        :param alpha: float (0-1], characterizes the grid size
        :param c: int (>0), determines number of grid cells
        :param est: str, estimator type (only 'mic_e' supported)
        """
        self.alpha = alpha
        self.c = c
        self.est = est
        self.x = None
        self.y = None
        self.mic = None
        self.tic = None

    def compute_score(self, x, y):
        #Core MIC computation between two vectors
        x = rankdata(x, method='average')
        y = rankdata(y, method='average')
        n = len(x)
        B = int(self.c * (n ** self.alpha))
        # Find the optimal grid partitioning and identify the max MIC.
        mic_val = 0
        tic_val = 0
        for b in range(2, B + 1):
            # Create grid bins
            x_bins = np.linspace(0, n, b + 1)[1:-1]
            y_bins = np.linspace(0, n, b + 1)[1:-1]

            # Calculate mutual information for this grid
            mi = self._mutual_information(x, y, x_bins, y_bins)
            norm = log(min(b, b))

            # Identify max MIC and compute TIC
            if mi / norm > mic_val:
                mic_val = mi / norm
            tic_val += mi / log(1 + b * b)

        self.mic = mic_val
        self.tic = tic_val / (B - 1)

    def _mutual_information(self, x, y, x_bins, y_bins):
        # Calculate mutual information between discretized variables
        # Discretize using bins
        x_disc = np.digitize(x, x_bins) 
        y_disc = np.digitize(y, y_bins) 

        # Joint probability matrix
        joint = np.histogram2d(x_disc, y_disc, bins=(len(x_bins) + 1, len(y_bins) + 1))[0]
        joint /= joint.sum()

        # Calculate MI using log2 for bits
        mi = 0.0
        for i in range(joint.shape[0]):      
            for j in range(joint.shape[1]):
                if joint[i, j] > 0:
                    mi += joint[i, j] * log(joint[i, j] /(joint[i, :].sum() * joint[:, j].sum()))
        return mi

    def mic(self):
        return self.mic

    def tic(self, norm=True):
        return self.tic if norm else self.tic * log(1 + B * B)

def parallel_process(X, i, j, alpha, c, est):
    """
    Args:
        i: Index of the first gene.
        j: Index of the second gene.
        df: The gene expression DataFrame.
    Returns:
        tuple: (i, j, correlation_value1, correlation_value2).
    """
    mine = MINE(alpha=alpha, c=c, est=est)
    mine.compute_score(X[i, :], X[j, :])
    return (i, j, mine.mic, mine.tic)

def par_pstats(X, alpha, c=1, est="mic"):
    """
    Args:
        X: The gene expression DataFrame(genes as rows, cells as columns).
    Returns:
        tuple: The matrix for value1, and value2 for calculations for different gene pairs.
    """
    n_genes = X.shape[0]
    # 1. Create a list of all combinations of i and j
    combinations = [(i, j) for i in range(n_genes) for j in range(i + 1, n_genes)]

    # 2. Initialize the results matrix
    mic_matrix = np.zeros((n_genes, n_genes))
    tic_matrix = np.zeros((n_genes, n_genes))

    # 3. Use multiprocessing to parallelize the loop
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        results = pool.starmap(parallel_process, [(X, i, j, 0.6, 1, "mic") for i, j in combinations])

    # 4. Collect the results and assemble the matrices
    for i, j, value1, value2 in results:
        mic_matrix[i, j] = value1
        mic_matrix[j, i] = value1
        tic_matrix[i, j] = value2
        tic_matrix[j, i] = value2
    return mic_matrix, tic_matrix

def pstats(X, alpha=0.6, c=1, est="mic_e"):
    # Pairwise MIC matrix for all gene pairs
    n_genes = X.shape[0]
    mic_matrix = np.zeros((n_genes, n_genes))
    tic_matrix = np.zeros((n_genes, n_genes))

    for i, j in combinations(range(n_genes), 2):
        mine = MINE(alpha=alpha, c=c, est=est)
        mine.compute_score(X[i,:], X[j,:])
        mic_matrix[i, j] = mic_matrix[j, i] = mine.mic
        tic_matrix[i, j] = tic_matrix[j, i] = mine.tic
    return mic_matrix, tic_matrix

def significant_mic_gmm(mic_matrix):
    """
    Identify significant MIC values using distribution-based thresholds
    Parameters:
        mic_matrix (pd.DataFrame): MIC values matrix
        method (str): 'gmm' (mixture model) or 'quantile'
        q (float): Quantile threshold (for 'quantile' method)
    Returns:
        pd.DataFrame: Boolean matrix indicating significant pairs
    """
    # Extract upper triangle values (excluding diagonal)
    upper_tri = mic_matrix.values[np.triu_indices_from(mic_matrix, k=1)]

    # Gaussian Mixture Model with 2 components
    gmm = GaussianMixture(n_components=2, random_state=0)
    gmm.fit(upper_tri.reshape(-1, 1))

    # Get threshold as mean + 3SD of background component
    means = gmm.means_.flatten()
    covs = np.sqrt(gmm.covariances_.flatten())
    bg_component = np.argmin(means)
    threshold = means[bg_component] + 3 * covs[bg_component]

    # Create significance matrix
    sig_matrix = mic_matrix.copy()
    #np.fill_diagonal(sig_matrix.values, False). This should be closed to allow setting dianogal=1.0
    sig_matrix = sig_matrix > threshold
    return sig_matrix


########################################################################################################################
# Auxiliary Functions - Read KEGG files and Wikipathway files
def read_kegg(geneFile, pathFile, linkFile):
    # Read genes and pathways in KEGG.
    gene_symbols = {}
    kegg_pathway = {}
    kegg_link    = {}
    gene_id_set  = set()

    with open(geneFile) as f:
        for line in f:
            cols = line.strip('\r\n').split('\t')
            gene_id = cols[0]
            descri = cols[-1]
            if (',' in descri) or (';' in descri):
                symbol_descri = descri.split(';')[0]
                symbols = symbol_descri.split(", ")
                if gene_id not in gene_symbols:
                    gene_symbols[gene_id] = set()
                for symbol in symbols:
                    gene_symbols[gene_id].add(symbol)
                gene_id_set.add(gene_id)

    with open(pathFile) as f:
        for line in f:
            pathID, pathwayName = line.strip('\r\n').split('\t')
            pathID = pathID.replace("path:", '')
            pathwayName = pathwayName[:pathwayName.rfind(' -')]
            kegg_pathway[pathID] = pathwayName

    with open(linkFile) as f:
        for line in f:
            pathID, gene_id = line.strip('\r\n').split('\t')
            pathID = pathID.replace("path:", '')
            if pathID in kegg_pathway:
                if gene_id in gene_id_set:
                    if pathID not in kegg_link:
                        kegg_link[pathID] = set()
                    kegg_link[pathID].add(gene_id)

    totalGeneNum = len(gene_id_set)
    kegg_gene = {}
    for gene_id in gene_symbols:
        for gene_symbol in gene_symbols[gene_id]:
            kegg_gene[gene_symbol] = gene_id

    return kegg_gene, kegg_pathway, kegg_link, totalGeneNum
	
def read_wikip(wikiDir, kegg_gene):
    # Read genes and pathways in Wikipathways.
    wiki_pathway = {}
    wiki_link = {}

    pathway_files = os.listdir(wikiDir)
    pathway_comp = re.compile(' Name="(.+?)" ')
    gene_comp = re.compile(' TextLabel="(.+?)" ')
    for pathway_file in pathway_files:
        pathway_id = pathway_file.split('_')[-2]
        with open(wikiDir + pathway_file, encoding='utf-8') as f:
            for line in f:
                pathway_re = re.search(pathway_comp, line)
                if pathway_re:
                    pathway = pathway_re.group(1)
                    wiki_pathway[pathway_id] = pathway_re.group(1)
                gene_re = re.search(gene_comp, line)
                if gene_re:
                    symbol = gene_re.group(1)
                    if symbol in kegg_gene:
                        gene_id = kegg_gene[symbol]
                        if pathway_id in wiki_pathway:
                            if pathway_id not in wiki_link:
                                 wiki_link[pathway_id] = set()
                            wiki_link[pathway_id].add(gene_id)

    return wiki_pathway, wiki_link

########################################################################################################################
# 1st Step - Computing lncRNA-lncRNA correlations
def compute_lnc_lnc_corr(df_gene_sample_exp):
    df_lncRNA_sample_exp = df_gene_sample_exp.copy(deep=True)
    df_lncRNA_sample_exp = df_lncRNA_sample_exp.drop(df_lncRNA_sample_exp[df_lncRNA_sample_exp.iloc[:,1] != 'lncRNA'].index)
    df_lncRNA_sample_exp.drop(["genetype"], axis=1, inplace=True)

    # drop duplicate genes(keep the one with higher zero_prop) and drop genes with zero_prop>0.9
    df = df_lncRNA_sample_exp.copy(deep=True)
    # Calculate zero_prop for each gene and remove genes with >90% zeros
    zero_proportions = (df.iloc[:, 1:] == 0).mean(axis=1)
    filtered_df = df[zero_proportions <= arg_zeroP].copy()
    # Calculate zero counts for remaining genes and remove duplicate genes
    zero_counts = (filtered_df.iloc[:, 1:] == 0).sum(axis=1)
    # Add "zero counts" as the (temporary) last column
    filtered_df['_temp_zero_count'] = zero_counts
    # Find indices of minimum zero count for each gene group
    min_indices = filtered_df.groupby(filtered_df.columns[0])['_temp_zero_count'].idxmin()
    # Get de_duplicated dataFrame and clean up
    deduped_df = filtered_df.loc[min_indices].drop(columns=['_temp_zero_count'])
    # Reset index if desired (the above removal may change index)
    deduped_df.reset_index(drop=True, inplace=True)

    # Get the remaining lncRNAs
    df_lncRNA_sample_exp = deduped_df.copy(deep=True)
    list_lncRNA = df_lncRNA_sample_exp.iloc[:,0]
    # Convert df to numpy array (and exclude gene IDs)
    df_lncRNA_sample_exp.set_index(['GENE'], inplace=True)
    #nd_lncRNA_sample_exp = df_lncRNA_sample_exp.values
    nd_lncRNA_sample_exp = df_lncRNA_sample_exp.values.astype(np.float32)

    # Computing correlation using either the normal or parallel version
    #nd_lnc_lnc_mic_matrix, nd_lnc_lnc_tic_matrix = pstats(nd_lncRNA_sample_exp, alpha=0.6, c=1, est="mic_e")
    nd_lnc_lnc_mic_matrix, nd_lnc_lnc_tic_matrix = par_pstats(nd_lncRNA_sample_exp, alpha=arg_alpha, c=1, est="mic_e")

    # Get significant correlations
    np.fill_diagonal(nd_lnc_lnc_mic_matrix, 1.0)
    np.fill_diagonal(nd_lnc_lnc_tic_matrix, 1.0)

    df_lnc_lnc_mic_matrix = pd.DataFrame(nd_lnc_lnc_mic_matrix)
    df_lnc_lnc_tic_matrix = pd.DataFrame(nd_lnc_lnc_tic_matrix)

    df_lnc_lnc_mic_matrix.index = list_lncRNA.tolist()
    df_lnc_lnc_tic_matrix.index = list_lncRNA.tolist()
    df_lnc_lnc_mic_matrix.columns = list_lncRNA.tolist()
    df_lnc_lnc_tic_matrix.columns = list_lncRNA.tolist()

    sig_lnc_lnc_mic_gmm = significant_mic_gmm(df_lnc_lnc_mic_matrix)
    sig_lnc_lnc_tic_gmm = significant_mic_gmm(df_lnc_lnc_tic_matrix)

    sig_lnc_lnc_corr  = sig_lnc_lnc_mic_gmm | sig_lnc_lnc_tic_gmm
    sig_lnc_lnc_corr2 = sig_lnc_lnc_mic_gmm & sig_lnc_lnc_tic_gmm

    # Save df_lncRNA_lncRNA_corr for manually checking the overlap of the two
    sig_lnc_lnc_corr.to_csv(arg_OutFold + '/sig_lnc_lnc_corr_OR.csv', index=True)
    sig_lnc_lnc_corr2.to_csv(arg_OutFold + '/sig_lnc_lnc_corr_AND.csv', index=True)

    return sig_lnc_lnc_corr, list_lncRNA

########################################################################################################################
# 2nd Step - Generate lncRNA sets (each lncRNA's correlated lncRNAs)
def generate_lncRNA_sets(list_lncRNA, sig_lnc_lnc_corr):
    dict_lncRNAsets = {}
    for row in range(0, len(list_lncRNA)):
        deleted = []
        for col in range(0, len(list_lncRNA)):
            if  sig_lnc_lnc_corr.iat[row, col] == False:
                deleted.append(list_lncRNA[col])

        for i in range(0, len(list_lncRNA)):
            if i != row and list_lncRNA[i] not in deleted:
                for j in range(0, len(list_lncRNA)):
                    if list_lncRNA[j] not in deleted:
                        if sig_lnc_lnc_corr.iat[i, j] == False:
                           deleted.append(list_lncRNA[j])

        lncRNAset = set(list_lncRNA) - set(deleted)
        if list_lncRNA[row] not in deleted:
            lncRNAlist = sorted(list(lncRNAset))
        else:
            lncRNAlist = []
        dict_lncRNAsets.update({list_lncRNA[row]: lncRNAlist})

    with open(arg_OutFold + '/dict_lncRNAsets.txt', 'w') as file_lncRNAsets:
        for key, value in dict_lncRNAsets.items():
            file_lncRNAsets.write(f"{key}: {value}\n")

    return dict_lncRNAsets

########################################################################################################################
# 3rd Step - Compute lncRNA-target correlations
def compute_lnc_target_corr(list_lncRNA, df_gene_sample_exp):
    df = df_gene_sample_exp.copy(deep=True)

    list_gene = df['GENE']
    list_type = df['genetype']

    targetlist = []
    llist_gene = list_gene.tolist()
    llist_type = list_type.tolist()
    llist_lncRNA = list_lncRNA.tolist()
    for j in range(len(llist_gene)):
        if llist_type[j] != 'TF' and llist_gene[j] not in llist_lncRNA and llist_gene[j] not in list_markers:
            targetlist.append(j)
    df.drop(index=targetlist, inplace=True)

    # Drop duplicate genes and genes with zero_prop>arg_zeroP
    # Calculate zero_prop for each gene (exclude 1st/GENE column) and remove genes with >90% zeros
    zero_proportions = (df.iloc[:, 1:] == 0).mean(axis=1)
    filtered_df = df[zero_proportions <= arg_zeroP].copy()
    # Calculate zero counts for remaining genes and remove duplicate genes
    zero_counts = (filtered_df.iloc[:, 1:] == 0).sum(axis=1)
    # Add "zero counts" as the (temporary) last column
    filtered_df['_temp_zero_count'] = zero_counts
    # Find indices of minimum zero count for each gene group
    min_indices = filtered_df.groupby(filtered_df.columns[0])['_temp_zero_count'].idxmin()
    # Get de_duplicated DataFrame and clean up
    deduped_df = filtered_df.loc[min_indices].drop(columns=['_temp_zero_count'])
    list_gene = deduped_df['GENE']
    list_type = deduped_df['genetype']

    deduped_df.set_index(['GENE'], inplace=True)
    deduped_df.drop(columns=['genetype'], inplace=True)
    nd_lnc_target_exp = deduped_df.values

    # Compute correlations using either the normal or parallel version
    #nd_lnc_target_mic_matrix, nd_lnc_target_tic_matrix = pstats(nd_lnc_target_exp, alpha=arg_alpha, c=1, est="mic_e")
    nd_lnc_target_mic_matrix, nd_lnc_target_tic_matrix = par_pstats(nd_lnc_target_exp, alpha=arg_alpha, c=1, est="mic_e")

    np.fill_diagonal(nd_lnc_target_mic_matrix, 1.0)
    np.fill_diagonal(nd_lnc_target_tic_matrix, 1.0)

    df_lnc_target_mic_matrix = pd.DataFrame(nd_lnc_target_mic_matrix)
    df_lnc_target_tic_matrix = pd.DataFrame(nd_lnc_target_tic_matrix)

    df_lnc_target_mic_matrix.index = list_gene.tolist()
    df_lnc_target_tic_matrix.index = list_gene.tolist()
    df_lnc_target_mic_matrix.columns = list_gene.tolist()
    df_lnc_target_tic_matrix.columns = list_gene.tolist()
    
    sig_lnc_target_mic_gmm = significant_mic_gmm(df_lnc_target_mic_matrix)
    sig_lnc_target_tic_gmm = significant_mic_gmm(df_lnc_target_tic_matrix)
 
    sig_lnc_target_corr  = sig_lnc_target_mic_gmm | sig_lnc_target_tic_gmm
    sig_lnc_target_corr2 = sig_lnc_target_mic_gmm & sig_lnc_target_tic_gmm

    targetlist=[]
    llist_gene = list_gene.tolist()
    llist_type = list_type.tolist()
    for j in range(len(llist_gene)):
        if llist_type[j] != 'lncRNA':
            targetlist.append(j)
    sig_lnc_target_corr.drop(sig_lnc_target_corr.columns[targetlist], axis=1, inplace=True)
    sig_lnc_target_corr2.drop(sig_lnc_target_corr2.columns[targetlist], axis=1, inplace=True)

    sig_lnc_target_corr.to_csv(arg_OutFold + '/sig_lnc_target_corr_OR.csv', index=True)
    sig_lnc_target_corr2.to_csv(arg_OutFold + '/sig_lnc_target_corr_AND.csv', index=True)

    return sig_lnc_target_corr, list_lncRNA

########################################################################################################################
# 4th Step - Generate lncRNA-target sets upon both df_lncRNA_gene_DBS and df_lncRNA_gene_corr
def generate_lncRNA_target(df_lncRNA_gene_DBS, df_lncRNA_gene_corr, dict_lncRNAsets):
    list_lncRNACORR = df_lncRNA_gene_corr.columns
    list_lncRNADBS = df_lncRNA_gene_DBS.columns

    dict_lncRNAtarget_DBS = {}
    dict_lncRNAtarget_CORR = {}
    dict_lncRNAtarget_shared = {}

    rows1 = df_lncRNA_gene_corr.shape[0]
    cols1 = df_lncRNA_gene_corr.shape[1]
    rows2 = df_lncRNA_gene_DBS.shape[0]
    cols2 = df_lncRNA_gene_DBS.shape[1]

    for key, value in dict_lncRNAsets.items():
        set_tmpCORR = set()
        set_tmpDBS = set()

        if key in list_lncRNACORR and value != []:
            for lncRNA in value:
                if lncRNA in df_lncRNA_gene_corr.columns:
                    for target in df_lncRNA_gene_corr.index:
                        if df_lncRNA_gene_corr[lncRNA][target] == True:
                            set_tmpCORR.add(target)
        dict_lncRNAtarget_CORR.update({key: set_tmpCORR})

        if key in list_lncRNADBS and value != []:
            for lncRNA in value:
                if lncRNA in df_lncRNA_gene_DBS.columns:
                    for target in df_lncRNA_gene_DBS.index:
                        if df_lncRNA_gene_DBS[lncRNA][target] > arg_LncRNADBSCutoff:
                            set_tmpDBS.add(target)
        dict_lncRNAtarget_DBS.update({key: set_tmpDBS})

        shared = set_tmpCORR.intersection(set_tmpDBS)
        dict_lncRNAtarget_shared.update({key: shared})

    with open(arg_OutFold + '/dict_lncRNAtarget_shared.txt', 'w') as file_lncRNAtarget_shared:
        for key, value in dict_lncRNAtarget_shared.items():
            file_lncRNAtarget_shared.write(f"{key}: {value}\n")

    return dict_lncRNAtarget_shared

########################################################################################################################
# 5th Step  - Compute TF-TF correlations
def compute_TF_TF_corr(df_gene_sample_exp):
    df_TF_sample_exp = df_gene_sample_exp.copy(deep=True)
    df_TF_sample_exp = df_TF_sample_exp.drop(df_TF_sample_exp[df_TF_sample_exp.iloc[:,1] != 'TF'].index)
    df_TF_sample_exp.drop(["genetype"], axis=1, inplace=True)

    # drop duplicate genes(keep the one with higher zero_prop) and drop genes with zero_prop>0.9
    df = df_TF_sample_exp.copy(deep=True)
    # Calculate zero_prop for each gene and remove genes with >90% zeros
    zero_proportions = (df.iloc[:, 1:] == 0).mean(axis=1)
    filtered_df = df[zero_proportions <= arg_zeroP].copy()
    # Calculate zero counts for remaining genes and remove duplicate genes
    zero_counts = (filtered_df.iloc[:, 1:] == 0).sum(axis=1)
    # Add "zero counts" as the (temporary) last column
    filtered_df['_temp_zero_count'] = zero_counts
    # Find indices of minimum zero count for each gene group
    min_indices = filtered_df.groupby(filtered_df.columns[0])['_temp_zero_count'].idxmin()
    # Get de_duplicated dataFrame and clean up
    deduped_df = filtered_df.loc[min_indices].drop(columns=['_temp_zero_count'])
    # Reset index if desired (the above removal may change index)
    deduped_df.reset_index(drop=True, inplace=True)

    # Get the remaining TFs
    df_TF_sample_exp = deduped_df.copy(deep=True)
    list_TF = df_TF_sample_exp.iloc[:,0]
    # Convert df to numpy array (and exclude gene IDs)
    df_TF_sample_exp.set_index(['GENE'], inplace=True)
    #nd_TF_sample_exp = df_TF_sample_exp.values
    nd_TF_sample_exp = df_TF_sample_exp.values.astype(np.float32)

    # Computing correlation using either the normal or parallel version
    #nd_TF_TF_mic_matrix, nd_TF_TF_tic_matrix = pstats(nd_TF_sample_exp, alpha=0.6, c=1, est="mic_e")
    nd_TF_TF_mic_matrix, nd_TF_TF_tic_matrix = par_pstats(nd_TF_sample_exp, alpha=arg_alpha, c=1, est="mic_e")

    # Get significant correlations
    np.fill_diagonal(nd_TF_TF_mic_matrix, 1.0)
    np.fill_diagonal(nd_TF_TF_tic_matrix, 1.0)

    df_TF_TF_mic_matrix = pd.DataFrame(nd_TF_TF_mic_matrix)
    df_TF_TF_tic_matrix = pd.DataFrame(nd_TF_TF_tic_matrix)

    df_TF_TF_mic_matrix.index = list_TF.tolist()
    df_TF_TF_tic_matrix.index = list_TF.tolist()
    df_TF_TF_mic_matrix.columns = list_TF.tolist()
    df_TF_TF_tic_matrix.columns = list_TF.tolist()

    sig_TF_TF_mic_gmm = significant_mic_gmm(df_TF_TF_mic_matrix)
    sig_TF_TF_tic_gmm = significant_mic_gmm(df_TF_TF_tic_matrix)

    sig_TF_TF_corr  = sig_TF_TF_mic_gmm | sig_TF_TF_tic_gmm
    sig_TF_TF_corr2 = sig_TF_TF_mic_gmm & sig_TF_TF_tic_gmm

    # Save df_TF_TF_corr for manually checking the overlap of the two
    sig_TF_TF_corr.to_csv(arg_OutFold + '/sig_TF_TF_corr_OR.csv', index=True)
    sig_TF_TF_corr2.to_csv(arg_OutFold + '/sig_TF_TF_corr_AND.csv', index=True)

    return sig_TF_TF_corr, list_TF

########################################################################################################################
# 6th Step  - Generate TF sets (each TF's correlated TFs)
def generate_TF_sets(list_TF, sig_TF_TF_corr):
    dict_TFsets = {}
    for row in range(0, len(list_TF)):
        deleted = []
        for col in range(0, len(list_TF)):
            if  sig_TF_TF_corr.iat[row, col] == False:
                deleted.append(list_TF[col])

        for i in range(0, len(list_TF)):
            if i != row and list_TF[i] not in deleted:
                for j in range(0, len(list_TF)):
                    if list_TF[j] not in deleted:
                        if sig_TF_TF_corr.iat[i, j] == False:
                           deleted.append(list_TF[j])

        TFset = set(list_TF) - set(deleted)
        if list_TF[row] not in deleted:
            TFlist = sorted(list(TFset))
        else:
            TFlist = []
        dict_TFsets.update({list_TF[row]: TFlist})

    with open(arg_OutFold + '/dict_TFsets.txt', 'w') as file_TFsets:
        for key, value in dict_TFsets.items():
            file_TFsets.write(f"{key}: {value}\n")

    return dict_TFsets

########################################################################################################################
# 7th Step  - Compute TF-target correlations
def compute_TF_target_corr(list_TF, df_gene_sample_exp):
    df = df_gene_sample_exp.copy(deep=True)

    list_gene = df['GENE']
    list_type = df['genetype']

    targetlist = []
    llist_gene = list_gene.tolist()
    llist_type = list_type.tolist()
    llist_TF   = list_TF.tolist()
    for j in range(len(llist_gene)):
        if llist_type[j] != 'lncRNA' and llist_gene[j] not in llist_TF and llist_gene[j] not in list_markers:
            targetlist.append(j)
    df.drop(index=targetlist, inplace=True)

    # Drop duplicate genes and genes with zero_prop>arg_zeroP
    # Calculate zero_prop for each gene (exclude 1st/GENE column) and remove genes with >90% zeros
    zero_proportions = (df.iloc[:, 1:] == 0).mean(axis=1)
    filtered_df = df[zero_proportions <= arg_zeroP].copy()
    # Calculate zero counts for remaining genes and remove duplicate genes
    zero_counts = (filtered_df.iloc[:, 1:] == 0).sum(axis=1)
    # Add "zero counts" as the (temporary) last column
    filtered_df['_temp_zero_count'] = zero_counts
    # Find indices of minimum zero count for each gene group
    min_indices = filtered_df.groupby(filtered_df.columns[0])['_temp_zero_count'].idxmin()
    # Get de_duplicated DataFrame and clean up
    deduped_df = filtered_df.loc[min_indices].drop(columns=['_temp_zero_count'])
    list_gene = deduped_df['GENE']
    list_type = deduped_df['genetype']

    deduped_df.set_index(['GENE'], inplace=True)
    deduped_df.drop(columns=['genetype'], inplace=True)
    nd_TF_target_exp = deduped_df.values

    # Compute correlations using either the normal or parallel version
    #nd_TF_target_mic_matrix, nd_TF_target_tic_matrix = pstats(nd_TF_target_exp, alpha=arg_alpha, c=1, est="mic_e")
    nd_TF_target_mic_matrix, nd_TF_target_tic_matrix = par_pstats(nd_TF_target_exp, alpha=arg_alpha, c=1, est="mic_e")

    np.fill_diagonal(nd_TF_target_mic_matrix, 1.0)
    np.fill_diagonal(nd_TF_target_tic_matrix, 1.0)

    df_TF_target_mic_matrix = pd.DataFrame(nd_TF_target_mic_matrix)
    df_TF_target_tic_matrix = pd.DataFrame(nd_TF_target_tic_matrix)

    df_TF_target_mic_matrix.index = list_gene.tolist()
    df_TF_target_tic_matrix.index = list_gene.tolist()
    df_TF_target_mic_matrix.columns = list_gene.tolist()
    df_TF_target_tic_matrix.columns = list_gene.tolist()

    sig_TF_target_mic_gmm = significant_mic_gmm(df_TF_target_mic_matrix)
    sig_TF_target_tic_gmm = significant_mic_gmm(df_TF_target_tic_matrix)

    sig_TF_target_corr = sig_TF_target_mic_gmm | sig_TF_target_tic_gmm
    sig_TF_target_corr2 = sig_TF_target_mic_gmm & sig_TF_target_tic_gmm

    targetlist = []
    llist_gene = list_gene.tolist()
    llist_type = list_type.tolist()
    for j in range(len(llist_gene)):
        if llist_type[j] != 'TF':
            targetlist.append(j)
    sig_TF_target_corr.drop(sig_TF_target_corr.columns[targetlist], axis=1, inplace=True)
    sig_TF_target_corr2.drop(sig_TF_target_corr2.columns[targetlist], axis=1, inplace=True)

    sig_TF_target_corr.to_csv(arg_OutFold + '/sig_TF_target_corr_OR.csv', index=True)
    sig_TF_target_corr2.to_csv(arg_OutFold + '/sig_TF_target_corr_AND.csv', index=True)

    return sig_TF_target_corr, list_TF

########################################################################################################################
# 8th Step  - Generate TF target sets upon both df_TF_gene_DBS and df_TF_gene_corr
def generate_TF_target(df_TF_gene_DBS, df_TF_gene_corr, dict_TFsets):
    list_TFCORR = df_TF_gene_corr.columns
    list_TFDBS = df_TF_gene_DBS.columns

    dict_TFtarget_DBS = {}
    dict_TFtarget_CORR = {}
    dict_TFtarget_shared = {}

    rows1 = df_TF_gene_corr.shape[0]
    cols1 = df_TF_gene_corr.shape[1]
    rows2 = df_TF_gene_DBS.shape[0]
    cols2 = df_TF_gene_DBS.shape[1]

    for key, value in dict_TFsets.items():
        set_tmpCORR = set()
        set_tmpDBS = set()

        if key in list_TFCORR and value != []:
            for TF in value:
                if TF in df_TF_gene_corr.columns:
                    for target in df_TF_gene_corr.index:
                        if df_TF_gene_corr[TF][target] == True:
                            set_tmpCORR.add(target)
        dict_TFtarget_CORR.update({key: set_tmpCORR})

        if key in list_TFDBS and value != []:
            for TF in value:
                if TF in df_TF_gene_DBS.columns:
                    for target in df_TF_gene_DBS.index:
                        if df_TF_gene_DBS[TF][target] > arg_TFDBSCutoff:
                            set_tmpDBS.add(target)
        dict_TFtarget_DBS.update({key: set_tmpDBS})

        shared = set_tmpCORR.intersection(set_tmpDBS)
        dict_TFtarget_shared.update({key: shared})

    with open(arg_OutFold + '/dict_TFtarget_shared.txt', 'w') as file_TFtarget_shared:
        for key, value in dict_TFtarget_shared.items():
            file_TFtarget_shared.write(f"{key}: {value}\n")

    return dict_TFtarget_shared

########################################################################################################################
# 9th Step - Identify sets and merge sets whose lncRNAs/TFs have the same targets
def collect_modules(dict_lncRNAtarget_shared, dict_TFtarget_shared, dict_lncRNAsets, dict_TFsets):
    '''
    dict_lncRNAsets            : returned by generate_lncRNA_sets().
	                             each key is a lncRNA and each value (a list) is the lncRNA's partners.
    dict_lncRNAtarget_shared   : returned by generate_lncRNA_target().
	                             each key is a lncRNA and each value (a set) is the lncRNA's targets.
	                             (only lncRNAs in the DBS file have targets).

    dict_TFsets                : returned by generate_TF_sets().
	                             each key is a TF and each value (a list) is the TF's partners.
    dict_lncRNAtarget_shared   : returned by generate_TF_target().
	                             each key is a TF and each value (a set) is the TF's targets.
	                             (only TFs in the DBS file have targets).

    dict_module_regulators     : generated here.
                                 each key is a set (list) of targets and each value is the union of the regulator sets.
    dict_module_regulator_sets : generated here.
                                 each key is a set (list) of targets but each value is the list of the regulator sets.
    '''
    dict_module_regulators = {}
    dict_module_regulator_sets = {}

    # Sort target genes in modules for determining whether two modules are the same
    allModules = [tuple(sorted(module)) for module in list(dict_lncRNAtarget_shared.values()) + list(dict_TFtarget_shared.values()) if len(module) > arg_ModuleSize]

    for module in allModules:
        dict_module_regulators[module] = set()
        dict_module_regulator_sets[module] = []

    for lncRNA in dict_lncRNAtarget_shared:
        module = tuple(sorted(dict_lncRNAtarget_shared[lncRNA]))
        if module in dict_module_regulators:
            dict_module_regulators[module].add(lncRNA + '(*)')    # '*' indicates lncRNA
            regulator_set = dict_lncRNAsets[lncRNA]
            if regulator_set and regulator_set not in dict_module_regulator_sets[module]:
                dict_module_regulator_sets[module].append(regulator_set)
        #else:  # open this print if we want to check these modules.
        #    print("The module of this lncRNA is not in 'dict_module_regulators' because its target set < arg_ModuleSize)", lncRNA, module)

    for TF in dict_TFtarget_shared:
        module = tuple(sorted(dict_TFtarget_shared[TF]))
        if module in dict_module_regulators:
            dict_module_regulators[module].add(TF + '(#)')       # '#' indicates TF
            regulator_set = dict_TFsets[TF]
            if regulator_set and regulator_set not in dict_module_regulator_sets[module]:
                dict_module_regulator_sets[module].append(regulator_set)
        #else:  # open this print if we want to check these modules.
        #    print("The module of this TF is not in 'dict_module_regulators' because its target set < arg_ModuleSize", TF, module)

    with open(arg_OutFold + '/dict_module_regulator_sets.txt', 'w') as file_module_regulator_sets:
        for key, value in dict_module_regulator_sets.items():
            file_module_regulator_sets.write(f"{key}: {value}\n")
    with open(arg_OutFold + '/dict_module_regulators.txt', 'w') as file_module_regulators:
        for key, value in dict_module_regulators.items():
            file_module_regulators.write(f"{key}: {value}\n")

    return dict_module_regulators, dict_module_regulator_sets

########################################################################################################################
# 10th Step - Pathway enrichment analysis for module genes
def pathway_analysis(dict_module_regulators, kegg_gene, kegg_pathway, kegg_link, totalGeneNum):
    dict_module_kegg = {}
    results = []
    pvalues = []
    
    for module in dict_module_regulators:
        module_geneIDs = {}
        for gene in module:
            gene = gene.strip('*#&_')
            if gene in kegg_gene:
                geneID = kegg_gene[gene]
                module_geneIDs[geneID] = gene
                
        moduleSize = len(module_geneIDs)
        
        for pathID in kegg_link:
            interSet = set(module_geneIDs.keys()) & kegg_link[pathID]
            if interSet:
                hitNum = len(interSet)
                thisPathGeneNum = len(kegg_link[pathID])
                
                # Indentify significantly enriched pathways using hypergeometric distribution test
                pVal = hypergeom.sf(hitNum-1, totalGeneNum, thisPathGeneNum, moduleSize)
                if pVal<0.05:
                    pathName = kegg_pathway[pathID]
                    hitGenes = set([module_geneIDs[geneID] for geneID in interSet])
                    results.append([module, pathID, pathName, hitGenes, pVal])
                    pvalues.append(pVal)
                
    # Benjamini-Hochberg method for p-value correction
    pValues = np.asarray(pvalues)
    by_descend = pValues.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(pValues)) / np.arange(len(pValues), 0, -1)
    fdr_values = np.minimum(1, np.minimum.accumulate(steps * pValues[by_descend]))
    fdrs = []
    for i in range(len(fdr_values[by_orig])):
        fdrs.append(fdr_values[by_orig][i])
        
    for i in range(len(fdrs)):
        results[i].append(fdrs[i])
        
    for module, pathID, pathName, hitGenes, pVal, fdr in results:
        if fdr < arg_fdrCutoff:
            if module not in dict_module_kegg:
                dict_module_kegg[module] = [(fdr, pathID, pathName, hitGenes)]
            else:
                dict_module_kegg[module].append((fdr, pathID, pathName, hitGenes))

    return dict_module_kegg

########################################################################################################################
# 11th Step - Merge modules that share a high percentage of targets
def remove_redundancy(curr_module_regulators, curr_module_regulator_sets):
    while True:
        keys   = list(curr_module_regulator_sets.keys())
        n = len(keys)
        if n < 2:
            break

        # Create N×N dataframes for both key and value similarity
        df_key_ratio   = pd.DataFrame(0.0, index=range(n), columns=range(n))
        df_value_ratio = pd.DataFrame(0.0, index=range(n), columns=range(n))

        # Calculate all pairwise similarities
        for i, j in combinations(range(n), 2):
            # Key similarity
            key_i, key_j = set(keys[i]), set(keys[j])
            shared_key = len(key_i & key_j)
            key_ratio = shared_key / min(len(key_i), len(key_j))
            df_key_ratio.at[i, j] = round(key_ratio,2)
            df_key_ratio.at[j, i] = round(key_ratio,2)

            # Value similarity (only calculate when key ratio is max)
            if key_ratio == 1.0:
                val_i = set(curr_module_regulator_sets[keys[i]][0])
                val_j = set(curr_module_regulator_sets[keys[j]][0])
                shared_val = len(val_i & val_j)
                value_ratio = shared_val / min(len(val_i), len(val_j)) if min(len(val_i), len(val_j)) > 0 else 0
                df_value_ratio.at[i, j] = round(value_ratio,2)
                df_value_ratio.at[j, i] = round(value_ratio,2)

        # Find maximum key ratio
        max_key_ratio = df_key_ratio.max().max()
        if max_key_ratio < arg_moduleOverlap:
            break

        # Get all pairs with maximum key ratio
        max_pairs = [(i, j) for i in range(n) for j in range(i + 1, n)
                     if df_key_ratio.at[i, j] == max_key_ratio]

        # If multiple pairs have ratio=1.0, choose the one with the highest value similarity
        if max_key_ratio == 1.0 and len(max_pairs) > 1:
            max_val_ratio = -1
            best_pair = max_pairs[0]
            for i, j in max_pairs:
                if df_value_ratio.at[i, j] > max_val_ratio:
                    max_val_ratio = df_value_ratio.at[i, j]
                    #print("max_val_ratio=", max_val_ratio)
                    best_pair = (i, j)
        else:
            best_pair = max_pairs[0]

        # Perform the merge
        i, j = best_pair
        keyi, keyj = keys[i], keys[j]
        new_key   = tuple(sorted(set(keyi) | set(keyj)))
        new_value = [list(set(curr_module_regulator_sets[keyi][0]) | set(curr_module_regulator_sets[keyj][0]))]
        new_val2  = set(curr_module_regulators[keyi]) | set(curr_module_regulators[keyj])

        # Remove old items and add merged item in both dict
        del curr_module_regulator_sets[keyi]
        del curr_module_regulator_sets[keyj]
        curr_module_regulator_sets[new_key] = new_value

        del curr_module_regulators[keyi]
        del curr_module_regulators[keyj]
        curr_module_regulators[new_key]     = new_val2

    with open(arg_OutFold + '/dict_module_regulators_noR.txt', 'w') as tmpfile1:
        for key, value in curr_module_regulators.items():
            tmpfile1.write(f"{key}: {value}\n")

    with open(arg_OutFold + '/dict_module_regulator_sets_noR.txt', 'w') as tmpfile2:
        for key, value in curr_module_regulator_sets.items():
            tmpfile2.write(f"{key}: {value}\n")

    return curr_module_regulators, curr_module_regulator_sets


########################################################################################################################
# 12th Step - Write results to files in the given directory
def write_modules(dict_module_regulators, dict_module_regulator_sets, dict_module_kegg, dict_module_wiki):
    with open(arg_OutFold + '/main.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['ModuleID',
                         'LncRNA(*)/TF(#)',
                         'Regulator set',
                         'Target gene',
                         'KEGG pathway, pathwayID, and FDR',
                         'Hit genes of enriched KEGG pathways',
                         'WikiPathway, pathwayID, and FDR',
                         'Hit genes of enriched Wikipathways'])
        moduleID = 0
        for module in dict_module_regulators:
            moduleID += 1
            if module in dict_module_kegg:
                this_module_kegg = []
                this_module_kegg_hitGenes = set()
                for fdr, pathID, pathName, hitGenes in sorted(dict_module_kegg[module]):
                    this_module_kegg.append(pathName + ' (' + pathID + ', ' + str('{0:3f}'.format(fdr)) + ')')
                    this_module_kegg_hitGenes = this_module_kegg_hitGenes | hitGenes
                pathway_kegg = '; '.join(this_module_kegg)
                hitGenes_kegg = ', '.join(this_module_kegg_hitGenes)
            else:
                pathway_kegg = hitGenes_kegg = ''

            if module in dict_module_wiki:
                this_module_wiki = []
                this_module_wiki_hitGenes = set()
                for fdr, pathID, pathName, hitGenes in sorted(dict_module_wiki[module]):
                    this_module_wiki.append(pathName + ' (' + pathID + ', ' + str('{0:3f}'.format(fdr)) + ')')
                    this_module_wiki_hitGenes = this_module_wiki_hitGenes | hitGenes
                pathway_wiki = '; '.join(this_module_wiki)
                hitGenes_wiki = ', '.join(this_module_wiki_hitGenes)
            else:
                pathway_wiki = hitGenes_wiki = ''

            regulators = ' | '.join(dict_module_regulators[module])
            regulatory_sets = ' | '.join(
                [str(regulatory_set).replace("'", '') for regulatory_set in dict_module_regulator_sets[module]])
            targetGenes = ' | '.join(module)

            writer.writerow(['ME' + str(moduleID),
                             regulators,
                             regulatory_sets,
                             targetGenes,
                             pathway_kegg,
                             hitGenes_kegg,
                             pathway_wiki,
                             hitGenes_wiki])


# Generate and write files for module visulization using cytoscape.
def generate_write_cytosc_file():
    f_re = open(arg_OutFold + '/moduleEdge', 'w')
    f_re.write('\t'.join(['lncRNA', 'module', 'weight']) + '\n')
    # allLncs = set()
    # allTFs = set()
    with open(arg_OutFold + '/main' + '.csv') as f_mod:
        reader = csv.reader(f_mod)
        first_row = next(reader)
        for row in reader:
            moduleID = row[0]
            regulators = row[1].split(' | ')

            for regulator in regulators:
                regulator_symbol = regulator.split('(')[0]
                f_re.write('\t'.join([regulator_symbol, moduleID, 'NA']) + '\n')
    f_re.close()


def write_lnc_target_relation(dict_lncRNAtarget_shared):
    with open(arg_OutFold + '/lnc_target_log.txt', 'w') as file_lnc_target_log:
        for key, value in dict_lncRNAtarget_shared.items():
            file_lnc_target_log.write(f"{key}: {value}\n")


########################################################################################################################
def usage():
    print("""Parameters:
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
        """)
    print("""Example: python eGRAMv3R1.py 
        -exp             data/M_SPG_2025.csv
        -lncDBS          data/Mmarker_MlncRNA_DBS.csv
        -tfDBS           data/Mmarker_MTF_DBS.csv
        -lncCutoff       100
        -tfCutoff        10
        -moduleSize      50
        -fdr             0.01
        -zeroP           0.85
        -merge           y
        -moduleOverlap   0.98
        -outFold         M_SG_out
        -species         2
        -alpha           0.0
        """)
    print()

#python  .\eGRAMv3R1.py --exp data/M_SPG_2025.csv --lncDBS data/Mmarker_MlncRNA_DBS.csv --lncCutoff 36 --moduleSize 50 --fdr 0.01 --zeroP 0.85 --merge y --moduleOverlap 0.98 --outFold M-SG-1 --species 2 --alpha 0.0

########################################################################################################################
########################################################################################################################
if __name__ == '__main__':
    param = sys.argv[1:]
    try:
        opts, args = getopt.getopt(param, '-h', ['exp=', 'lncDBS=', 'lncCutoff=', 'moduleSize=', 'fdr=', 'zeroP=', 'merge=', 'moduleOverlap=', 'outFold=', 'species=', 'alpha='])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit(2)
        elif opt == '--exp':
            df_gene_sample_exp = pd.read_csv(arg)
        elif opt == '--lncDBS':
            df_lncRNA_gene_DBS = pd.read_csv(arg, index_col=[0])
        elif opt == '--tfDBS':
            df_TF_gene_DBS = pd.read_csv(arg, index_col=[0])
        elif opt == '--lncCutoff':
            arg_LncRNADBSCutoff = int(arg)
        elif opt == '--tfCutoff':
            arg_TFDBSCutoff = int(arg)
        elif opt == '--moduleSize':
            arg_ModuleSize = int(arg)
        elif opt == '--fdr':
            arg_fdrCutoff = float(arg)
        elif opt == '--zeroP':
            arg_zeroP = float(arg)
        elif opt == '--merge':
            arg_RemoveRedund = arg
        elif opt == '--moduleOverlap':
            arg_moduleOverlap = float(arg)
        elif opt == '--outFold':
            arg_OutFold = arg
        elif opt == '--species':
            arg_SpeciesPath = int(arg)
        elif opt == '--alpha':
            arg_alpha = float(arg)

    if 'tfDBS' not in dir():
        df_TF_gene_DBS = pd.DataFrame()

    if 'arg_LncRNADBSCutoff' not in dir():
        arg_LncRNADBSCutoff = 36
    if 'arg_TFDBSCutoff' not in dir():
        arg_TFDBSCutoff = 10
    if 'arg_ModuleSize' not in dir():
        arg_ModuleSize = 2                 #50
    if 'arg_fdrCutoff' not in dir():
        arg_fdrCutoff = 0.01
    if 'arg_zeroP' not in dir():
        arg_zeroP = 0.85
    if 'arg_RemoveRedund' not in dir():
        arg_RemoveRedund = 'y'
    if 'arg_moduleOverlap' not in dir():
        arg_moduleOverlap = 0.98
    if 'arg_OutFold' not in dir():
        arg_OutFold = "Output/"
    if 'arg_SpeciesPath' not in dir():
        arg_SpeciesPath = 2
    if 'arg_alpha' not in dir():
        arg_alpha = 0.0

########################################################################################################################
    if not os.path.exists(arg_OutFold):
        os.mkdir(arg_OutFold)

    if arg_SpeciesPath == 1:
        keggGeneFile = "PathwayAnnotation/keggGene-hsa"  
        keggPathFile = "PathwayAnnotation/keggPathway-hsa"  
        keggLinkFile = "PathwayAnnotation/keggLink-hsa"  
        wikiPath = "PathwayAnnotation/wikipathways-20240310-gpml-Homo_sapiens/"  
    elif arg_SpeciesPath == 2:
        keggGeneFile = "PathwayAnnotation/keggGene-mmu"
        keggPathFile = "PathwayAnnotation/keggPathway-mmu"
        keggLinkFile = "PathwayAnnotation/keggLink-mmu"
        wikiPath = 'PathwayAnnotation/wikipathways-20240310-gpml-Mus_musculus/'
    else:
        sys.exit(2)

    if arg_alpha == 0.0:   # Control time consumption of the MIC algorithm if sample size is large
        samplesize = df_gene_sample_exp.shape[1]
        if samplesize<500:
            arg_alpha=0.6
        elif samplesize<1000:
            arg_alpha=0.55
        elif samplesize<1500:
            arg_alpha=0.5
        else:
            arg_alpha=0.45

    if df_TF_gene_DBS.empty:
        TFflag = 0
    else:
        TFflag = 1

    kegg_gene, kegg_pathway, kegg_link, totalGeneNum = read_kegg(keggGeneFile, keggPathFile, keggLinkFile)
    wiki_pathway, wiki_link = read_wikip(wikiPath, kegg_gene)

########################################################################################################################
    # 1st Step - compute lncRNA-lncRNA correlations
    df_lncRNA_lncRNA_corr, list_lncRNA = compute_lnc_lnc_corr(df_gene_sample_exp)

    # 2nd Step - generate_lncRNA_sets
    dict_lncRNAsets = generate_lncRNA_sets(list_lncRNA, df_lncRNA_lncRNA_corr)

    # 3rd Step - compute lncRNA-gene correlations
    df_lncRNA_gene_corr, list_lncRNA   = compute_lnc_target_corr(list_lncRNA, df_gene_sample_exp)

    # 4th Step - generate shared lncRNA-target sets upon df_lncRNA_gene_DBS and df_lncRNA_gene_corr
    dict_lncRNAtarget_shared = generate_lncRNA_target(df_lncRNA_gene_DBS, df_lncRNA_gene_corr, dict_lncRNAsets)

########################################################################################################################
    if TFflag == 1:
        # 5th Step - compute TF-TF correlations
        df_TF_TF_corr, list_TF = compute_TF_TF_corr(df_gene_sample_exp)

        # 6th Step - generate_TF_sets
        dict_TFsets = generate_TF_sets(list_TF, sig_TF_TF_corr)

        # 7th Step - compute TF-gene correlations
        df_TF_gene_corr, list_TF = compute_TF_target_corr(list_TF, df_gene_sample_exp)

        # 8th Step - generate shared TF target sets upon df_TF_gene_DBS and df_TF_gene_corr
        dict_TFtarget_shared = generate_TF_target(df_TF_gene_DBS, df_TF_gene_corr, dict_TFsets)
    else:
        #dict_TFtarget_shared = {}
        #dict_TFsets = {}
        dict_TFtarget_shared = {('_Chd6', '_Celf2', '_Dock3'): []}
        dict_TFsets          = {'_Lpp': ['_Chd6', '_Celf2', '_Dock3']}


########################################################################################################################
    # 9th Step - identify modules and merge lncRNAs/TFs that regulate the same modules
    #if TFflag == 1:
    dict_module_regulators, dict_module_regulator_sets = collect_modules(dict_lncRNAtarget_shared, dict_TFtarget_shared, dict_lncRNAsets, dict_TFsets)
    #else:
    #    dict_module_regulators, dict_module_regulator_sets = collect_modules(dict_lncRNAtarget_shared, dict_lncRNAsets)

    # 10th Step - pathway enrichment analysis for module genes in df_path_module (KEGG gene annotation is used)
    dict_module_kegg = pathway_analysis(dict_module_regulators, kegg_gene, kegg_pathway, kegg_link, totalGeneNum)
    dict_module_wiki = pathway_analysis(dict_module_regulators, kegg_gene, wiki_pathway, wiki_link, totalGeneNum)

    # 11th Step - merge modules
    if arg_RemoveRedund == 'y':
        dict_module_regulators_noR, dict_module_regulator_sets_noR = remove_redundancy(dict_module_regulators, dict_module_regulator_sets)

    # 12th Step - write results to files
    write_modules(dict_module_regulators, dict_module_regulator_sets, dict_module_kegg, dict_module_wiki)
    write_lnc_target_relation(dict_lncRNAtarget_shared)
    generate_write_cytosc_file()
