#Import dependencies
from bs4 import BeautifulSoup
import requests as re
import json

def merge(control_samples, treated_samples): 
    res = {**control_samples, **treated_samples} 
    return res 

url = 'https://amp.pharm.mssm.edu/Harmonizome/api/1.0/gene_set/autophagy/GO+Biological+Process+Annotations'
page1 = re.get(url)
soup = str(BeautifulSoup(page1.text, 'html.parser'))
data = json.loads(soup)

list1 = data['associations']

gene_list = []
for i in range(0, (len(list1))):
    dict1 = list1[i-1]
    dict2 = list1[i]
    dict3 = merge(dict1, dict2)
    gene_list.append(dict3['gene']['symbol'])
    
gene_dict = {'GOBioProcAnnot': gene_list}