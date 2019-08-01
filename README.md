# Microarray Analysis Pipeline in Python

<b> Abstract: </b> Autophagy is a catabolic pathway in which cytoplasmic vesicles engulf long-lived cellular components to transport them to the lysosome for digestion. Autophagy is a critical process for maintaining cell homeostasis, and stress-inducing conditions inhibit or activate autophagy, which can lead to a wide range of pathologies. Utilizing the NIH Gene Expression Omnibus (GEO) database to obtain gene expression data from autophagy-related studies, I have constructed an analysis pipeline in Jupyter Notebook to determine differentially expressed genes (DEGs). With the user's input of the GEO accession ID and control and treated samples, or a file with relevant metadata, the pipeline performs normalization, dimensionality reduction, and the Characteristic Direction method to generate lists of upregulated and downregulated genes. The gene sets from over 100 studies were compiled with sets from databases including Gene Ontology and KEGG pathways; studies with significant overlap with these resources were utilized to determine a more extensive genetic signature for autophagy. Once the consensus signature was generated, the data was imported into the L1000FWD web application to establish novel drugs that could induce autophagy. Since autophagy is linked to conditions such as cancer and neurodegeneration, targeting this pathway could transform established treatments for disease.

### Documentation:

Data:
- RNA-seq data extracted from BioJupies (3)
- Rapamycin gene signature from CREEDS
- Autophagy gene signature predictions from Geneshot (6)
- example_metadata.csv: Example of metadata, focused on autophagy
- probe2gene.txt.zip: Text file used to map microarray probes to gene symbols

Notebooks:
- Analysis from CSV to Clustergrammer Matrix.ipynb: Takes CSV of GEO metadata, performs differential expression analysis, and exports as Jaccard similarity matrix (for Clustergrammer heatmap visualization)
- Analysis from CSV to GMT.ipynb: Takes CSV of GEO metadata, performs differential expression analysis and exports as GMT file
- Microarray Analysis Single Study.ipynb: Takes user input of accession id and control and treated samples to perform differential and enrichment analysis on gene expression data

Scripts:
- microarray_analysis_function.py: Function for performing microarray analysis, returns lists of up- and down-regulated genes
