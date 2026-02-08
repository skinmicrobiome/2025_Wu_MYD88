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

###FigS5###
df = pd.read_csv('metatble.csv', dtype={'subject_id': str})
df['mouseID'] = df['mouseID'] + '.' + df['timepoint']
df['mouseID'] = df['mouseID'].str.replace("ay ", '')
df['timepoint'] = df['timepoint'].str.replace("ay ", '')

group_info = (
    df[['mouseID', 'index', 'treatment_group']]
    .drop_duplicates()
    .set_index('mouseID')
)

group_dict = {
    'ft/ft': ['ft/ft.f1.D0', 'ft/ft.f2.D0', 'ft/ft.f3.D0', 'ft/ft.f4.D0'],
    'ft/ftXMyD88-/-': ['neo.f1.D0', 'neo.f2.D0', 'neo.f3.D0', 'neo.f4.D0', 'neo.f5.D0'],
    'ft/ftXMyD88-/- veh': ['veh.f1.D7', 'veh.f2.D7', 'veh.f3.D7', 'veh.f4.D7', 'veh.f5.D7'],
    'ft/ftXMyD88-/- neo': ['neo.f1.D7', 'neo.f2.D7', 'neo.f4.D7', 'neo.f5.D7']
}

lookup = {v: k for k, values in group_dict.items() for v in values}
group_info["group"] = group_info.index.map(lookup)

df1 = df.copy()
df1['rel_abun'] = (
    df1['new_est_reads'] /
    df1.groupby('mouseID')['new_est_reads'].transform('sum') * 100
)

df_pivot = df1.pivot(index='mouseID', columns='ASV', values='rel_abun')

print(df_pivot.sum(axis=1))

df_pivot = pd.merge(df_pivot, group_info, left_index=True, right_index=True).fillna(0)
df_pivot.to_csv('figs5.csv')