import pandas as pd
import numpy as np

# Read & filter to bacteria
df_otu = pd.read_csv('metatble.csv')
df_otu = df_otu[df_otu['kindom'] == 'Bacteria']

# metadata
info = df_otu[['treatment_group', 'timepoint', 'mouseID', 'index']].drop_duplicates()

# choose labels
phyla   = ['Firmicutes', 'Actinobacteriota', 'Bacteroidota', 
           'Proteobacteria', 'Verrucomicrobiota']
genera  = ['Streptococcus', 'Staphylococcus', 'Corynebacterium',
           'Bacteroides', 'Alistipes', 'Akkermansia']

df_otu['label']="Bacteria (Other)"

df_otu['label'] = np.where(df_otu['phylum'].isin([
'Firmicutes',
'Actinobacteriota',
'Bacteroidota',
'Proteobacteria', 
'Verrucomicrobiota']), df_otu['phylum'], df_otu['label'])

df_otu['label'] = np.where(df_otu['genus'].isin([
'Streptococcus', 
'Staphylococcus', 
'Corynebacterium', 
'Bacteroides', 
'Alistipes',
'Akkermansia']), df_otu['genus'], df_otu['label'])

update={"Actinobacteriota":"Actinomycetota",
       "Firmicutes":"Bacillota",
       "Proteobacteria":"Pseudomonadota"}
df_otu["label"] = df_otu["label"].replace(update)   

# relative abundance
df_otu['rel_abun'] =df_otu['new_est_reads']/ df_otu.groupby(['mouseID', 'timepoint'])['new_est_reads'].transform('sum')* 100
df = df_otu.groupby(['mouseID', 'timepoint', 'label'], as_index=False)['rel_abun'].sum()
df = pd.merge(df,info)

# generate
df.to_csv('fig4bd.csv', index=False)

#check and reads deposition
print('index:\n',set(df_otu['index'].to_list()))
print('label:\n',df_otu.label.unique())
print('mouseID:\n',df_otu['mouseID'].unique())
print('relative abundance:\n',df_otu.groupby(['mouseID','timepoint'])['rel_abun'].sum())
print('BioProject SampleId:\n',df_otu['SampleID'].unique())
