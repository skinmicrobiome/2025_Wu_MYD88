import pandas as pd
import numpy as np

phyla  = ['Firmicutes', 'Actinobacteriota', 'Bacteroidota',
          'Proteobacteria', 'Verrucomicrobiota']
genera = ['Streptococcus', 'Staphylococcus', 'Corynebacterium',
          'Bacteroides', 'Alistipes', 'Akkermansia']
rename = {'Actinobacteriota': 'Actinomycetota',
          'Firmicutes': 'Bacillota',
          'Proteobacteria': 'Pseudomonadota'}
group_dict = {
    'ft/ft':     ['ft/ft.f1.D0', 'ft/ft.f2.D0', 'ft/ft.f3.D0', 'ft/ft.f4.D0'],
    'neo.Day 0': ['neo.f1.D0',   'neo.f2.D0',   'neo.f3.D0',   'neo.f4.D0',   'neo.f5.D0'],
    'veh.Day 7': ['veh.f1.D7',   'veh.f2.D7',   'veh.f3.D7',   'veh.f4.D7',   'veh.f5.D7'],
    'neo.Day 7': ['neo.f1.D7',   'neo.f2.D7',   'neo.f4.D7',   'neo.f5.D7'],
}

df = pd.read_csv('metatble.csv')
df = df[df['kindom'] == 'Bacteria']

info = df[['treatment_group', 'timepoint', 'mouseID', 'index']].drop_duplicates()

# Fig5bd
df['label'] = 'Bacteria (Other)'
df['label'] = np.where(df['phylum'].isin(phyla), df['phylum'], df['label'])
df['label'] = np.where(df['genus'].isin(genera), df['genus'], df['label'])
df['label'] = df['label'].replace(rename) 
df['rel_abun'] = (df['new_est_reads']
                  / df.groupby(['mouseID', 'timepoint'])['new_est_reads'].transform('sum') * 100)

fig5bd = (df.groupby(['mouseID', 'timepoint', 'label'], as_index=False)['rel_abun']
          .sum()
          .merge(info))
fig5bd.to_csv('fig5bd.csv', index=False)

# FigS2 pivot
s2 = df.copy()
s2['mouseID']   = (s2['mouseID'] + '.' + s2['timepoint']).str.replace('ay ', '', regex=False)
s2['timepoint'] = s2['timepoint'].str.replace('ay ', '', regex=False)

lookup = {v: k for k, vs in group_dict.items() for v in vs}
s2_info = (s2[['mouseID', 'index', 'treatment_group']]
           .drop_duplicates()
           .set_index('mouseID')
           .assign(group=lambda x: x.index.map(lookup)))

s2['rel_abun'] = s2['new_est_reads'] / s2.groupby('mouseID')['new_est_reads'].transform('sum') * 100
s2_pivot = (s2.pivot(index='mouseID', columns='ASV', values='rel_abun')
              .merge(s2_info, left_index=True, right_index=True)
              .fillna(0))
s2_pivot.to_csv('figs2.csv')

# info
checks = {
    'index':                   sorted(df['index'].unique()),
    'labels':                  df['label'].unique(),
    'mouseIDs':                df['mouseID'].unique(),
    'rel_abun per mouse/time': df.groupby(['mouseID', 'timepoint'])['rel_abun'].sum(),
    'BioProject SampleIDs':    df['SampleID'].unique(),
}
for k, v in checks.items():
    print(f'\n{k}:\n{v}')
    