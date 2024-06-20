import pandas as pd
import numpy as np

df1 = pd.read_csv('uudet.csv')
df2 = pd.read_csv('vanhat.csv')

unique_to_df1 = df1[~df1['NiftiPath'].isin(df2['NiftiPath'])]
unique_to_df1.to_csv('uudet-uniikit.csv', index=False)

print('Koko:', unique_to_df1.shape[0])
