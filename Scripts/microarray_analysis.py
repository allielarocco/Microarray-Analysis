#Import required functions
import GEOparse
import pandas as pd
import numpy as np
import os
import json
from sklearn.preprocessing import quantile_transform
from sklearn.decomposition import PCA
import warnings
from scipy.stats import chi2
from scipy.stats.mstats import zscore
import time

#Time sleep to prevent crashes
time.sleep(1)
#Change this to your working directory
os.chdir('../Data')


#Define merge function for combining sample lists
def merge(dict1, dict2): 
        res = {**dict1, **dict2} 
        return res

#Define probe_dict
def probe_dict(probe2gene_file):
        dict1 = {}
        with open(probe2gene_file, 'r') as f:
            for line in f:
                line = line.strip()
                (platform, probe, symbol) = line.split()
                dict1[probe] = symbol
        return dict1

PROBE2GENE = probe_dict('../Data/probe2gene.txt')

#Define characteristic direction function
def chdir(data, sampleclass, genes, gamma=1., sort=True, calculate_sig=False, nnull=10, sig_only=False, norm_vector=True):
        data.astype(float)
        #sampleclass = np.array(map(int, sampleclass))
        m_non0 = sampleclass != 0
        m1 = sampleclass[m_non0] == 1
        m2 = sampleclass[m_non0] == 2

        data = data[:, m_non0]
        data = zscore(data)

        n1 = m1.sum() # number of controls
        n2 = m2.sum() # number of experiments

        meanvec = data[:,m2].mean(axis=1) - data[:,m1].mean(axis=1) 

        pca = PCA(n_components=None)
        pca.fit(np.array(data.T))

        cumsum = pca.explained_variance_ratio_
        keepPC = len(cumsum[cumsum > 0.001])

        v = pca.components_[0:keepPC].T
        r = pca.transform(data.T)[:,0:keepPC]
        dd = ( np.dot(r[m1].T,r[m1]) + np.dot(r[m2].T,r[m2]) ) / float(n1+n2-2)
        sigma = np.mean(np.diag(dd))

        shrunkMats = np.linalg.inv(gamma*dd + sigma*(1-gamma)*np.eye(keepPC))
        b = np.dot(v, np.dot(np.dot(v.T, meanvec), shrunkMats))

        if norm_vector:
            b /= np.linalg.norm(b)

        grouped = zip([abs(item) for item in b],b,genes)
        if sort:
            grouped = sorted(grouped,key=lambda x: x[0], reverse=True)

        if not calculate_sig: # return sorted b and genes.
            res = [(item[1],item[2]) for item in grouped]
            return res
        else: # generate a null distribution of chdirs
            nu = n1 + n2 - 2
            y1 = np.random.multivariate_normal(np.zeros(keepPC), dd, nnull).T * np.sqrt(nu / chi2.rvs(nu,size=nnull))
            y2 = np.random.multivariate_normal(np.zeros(keepPC), dd, nnull).T * np.sqrt(nu / chi2.rvs(nu,size=nnull))
            y = y2 - y1

            nullchdirs = []
            for col in y.T:
                bn = np.dot(np.dot(np.dot(v,shrunkMats), v.T), np.dot(col,v.T))
                bn /= np.linalg.norm(bn)
                bn = bn ** 2
                bn.sort()
                bn = bn[::-1] ## sort in decending order
                nullchdirs.append(bn)

            nullchdirs = np.array(nullchdirs).T
            nullchdirs = nullchdirs.mean(axis=1)
            b_s = b ** 2 
            b_s.sort()
            b_s = b_s[::-1] # sorted b in decending order
            relerr = b_s / nullchdirs ## relative error
            # ratio_to_null
            ratios = np.cumsum(relerr)/np.sum(relerr)- np.linspace(1./len(meanvec),1,len(meanvec))
            res = [(item[1],item[2], ratio) for item, ratio in zip(grouped, ratios)] 
            print('Number of significant genes: %s'%(np.argmax(ratios)+1))
            if sig_only:
                return res[0:np.argmax(ratios)+1]
            else:
                return res


#Microarray Analysis Pipeline
def micro_analysis(accession_id, control_samples, treated_samples):
    #Creating a dictionary of assigned control and treated samples 
      
    control_samples = { i : 'control' for i in control_samples }
    treated_samples = { i : 'treated' for i in treated_samples }
    all_samples = merge(control_samples, treated_samples)
    
    #Parse the GEO data using the Accession ID
    gse = GEOparse.get_GEO(geo=accession_id, destdir="./")
    #Create a list of samples to use in the development of the expression matrix
    list_samples = list(all_samples.keys())
    
    #Visualization of expression matrix
    pivoted_samples = gse.pivot_samples('VALUE')[list_samples]
    pivoted_samples.head()
    #Determine the total amount of probes used in the study
    pivoted_samples_average = pivoted_samples.median(axis=1)
    #Filtering out unexpressed probes
    expression_threshold = pivoted_samples_average.quantile(0.3)
    expressed_probes = pivoted_samples_average[pivoted_samples_average >= expression_threshold].index.tolist()
    
    #Redefine expression data using only the expressed probes
    exprsdata = gse.pivot_samples("VALUE").loc[expressed_probes]
    exprsdata = exprsdata.T
    #Deletes additional samples that aren't being analyzed
    exprsdata = exprsdata[exprsdata.index.isin(list_samples)]
    #Drop any probe columns where expression data is missing or negative
    exprsdata.dropna(axis = 1)
    
    #Quantile normalization of data
    rank_mean = exprsdata.stack().groupby(exprsdata.rank(method='first').stack().astype(int)).mean()
    exprsdata.rank(method='min').stack().astype(int).map(rank_mean).unstack().dropna(axis=1)
    #Making Dataframe of samples
    samplesDf = pd.DataFrame.from_dict(all_samples, orient = 'index', columns = ['type'])
    samplesDf.reset_index(inplace=True)
    
    #Transpose data matrix for sorting, index correlated to probe IDs
    exprsdata = exprsdata.T
    #Upload annotation file as dictionary

    #Reset index and replace with gene symbols, view as dataframe
    exprsdata = pd.DataFrame(exprsdata)
    exprsdata.index = exprsdata.index.astype(str, copy=False)
    exprsdata['symbol'] = exprsdata.index.to_series().map(PROBE2GENE)
    exprsdata.reset_index(inplace=True)
    data = exprsdata.set_index('symbol')
    
    #Drop probe id column
    data = data.drop('ID_REF', axis=1)
    #Drop rows that aren't associated with a particular gene symbol
    data = data.reset_index().dropna().set_index('symbol')
    
    #Utilize warning statements
    warnings.filterwarnings("ignore", category=DeprecationWarning) 
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    
    #Make sample classes, ensure that there is a distinction between control/treated samples
    data_cd = {}
    sample_classes = {}
    sample_class = np.zeros(data.shape[1], dtype=np.int32)
    sample_class[samplesDf['type'].values == 'control'] = 1
    sample_class[samplesDf['type'].values == 'treated'] = 2
    sample_classes = sample_class
    
    #CD results
    cd_res = chdir(data.values, sample_classes, data.index, gamma=.5, sort=False, calculate_sig=False)
    cd_coefs = np.array(list(map(lambda x: x[0], cd_res)))
    srt_idx = np.abs(cd_coefs).argsort()[::-1]
    cd_coefs = cd_coefs[srt_idx][:600]
    sorted_DEGs = data.index[srt_idx][:600]
    up_genes = dict(zip(sorted_DEGs[cd_coefs > 0], cd_coefs[cd_coefs > 0]))
    dn_genes = dict(zip(sorted_DEGs[cd_coefs < 0], cd_coefs[cd_coefs < 0]))
    data_cd['up'] = up_genes
    data_cd['dn'] = dn_genes
    
    #Retrieve up and down gene sets
    up_list = list(up_genes.keys())
    dn_list = list(dn_genes.keys())
    #Up genes and down genes
    return up_list, dn_list
    