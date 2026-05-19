import pandas as pd
import numpy as np

df = pd.read_excel('dataset/sod-1_data.xlsx')
print('Exact columns:', list(df.columns))
print()

print('Chemical patterns:')
for cp in sorted(df['Chemical_Pattern'].unique()):
    cnt = (df['Chemical_Pattern'] == cp).sum()
    print(f'  [{cp}] -> count: {cnt}')
print()

print('Modification values:')
for m in df['Modification'].unique():
    print(f'  [{m}]')
print()

print('Sample sequences (first 10):')
for i, s in enumerate(df['Sequence'].head(10)):
    print(f'  {i}: {s} (len={len(s)})')
print()

print('Inhibition distribution:')
print(df['Inhibition(%)'].describe())
print()

print('Sequence+ChemPattern combos:', df.groupby(['Sequence', 'Chemical_Pattern']).ngroups)
print('Seq+ChemPat+Linkage combos:', df.groupby(['Sequence', 'Chemical_Pattern', 'Linkage']).ngroups)
print()

# Check for duplicated rows
print('Fully duplicated rows:', df.duplicated().sum())
print('Duplicated sequences:', df['Sequence'].duplicated().sum())
print()

# Sample Location column
print('Sample Location values:')
for i, loc in enumerate(df['Location'].head(5)):
    print(f'  {i}: {loc}')
print()

# Sample Linkage_Location column
print('Sample Linkage_Location values:')
for i, ll in enumerate(df['Linkage_Location'].head(5)):
    print(f'  {i}: {ll}')
print()

# Primer_probe_set
print('Unique primer_probe_set:', df['Primer_probe_set'].nunique())
for pp in df['Primer_probe_set'].unique():
    cnt = (df['Primer_probe_set'] == pp).sum()
    print(f'  [{pp}] -> {cnt}')
print()

# Check sequence characters
from collections import Counter
all_chars = Counter()
for s in df['Sequence']:
    all_chars.update(s.upper())
print('Characters in sequences:', dict(all_chars))

# Check if seq_length matches actual length
df['computed_len'] = df['Sequence'].str.len()
mismatches = (df['seq_length'] != df['computed_len']).sum()
print(f'\nSeq length mismatches: {mismatches}')
if mismatches > 0:
    print(df[df['seq_length'] != df['computed_len']][['Sequence', 'seq_length', 'computed_len']].head(10))

# Density distribution
print('\nDensity distribution:')
print(df['Density(cells/well)'].value_counts())

# ASO_volume distribution
print('\nASO_volume distribution:')
print(df['ASO_volume(nM)'].value_counts().sort_index())

# Treatment_Period distribution
print('\nTreatment_Period distribution:')
print(df['Treatment_Period(hours)'].value_counts())

# Inhibition range check
print(f'\nInhibition min: {df["Inhibition(%)"].min()}, max: {df["Inhibition(%)"].max()}')
out_of_range = ((df['Inhibition(%)'] < 0) | (df['Inhibition(%)'] > 100)).sum()
print(f'Inhibition out of [0,100]: {out_of_range}')
