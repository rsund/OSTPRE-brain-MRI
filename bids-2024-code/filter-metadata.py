import pandas as pd

def contains_any_substring(cell, substrings):
    if isinstance(cell, str):
      return any(substring in cell for substring in substrings)
    return False

def get_interesting_rows(file_name):
  df = pd.DataFrame(pd.read_csv('{}.csv'.format(file_name)))

  wanted_words = "mpr|mprage|t1_se|sag|tra|T1|TFE"
  unwanted_words = "gd|GD|Gd|gD"
  column_to_check = "ProtocolName"
  columns_to_check = ['ProtocolName', 'SeriesDescription', 'PulseSequenceName']

  substrings = ['mpr', 'mprage', 't1_se', 'sag', 'tra', 'T1', 'TFE']
  unwanted_substrings = "gd|GD|Gd|gD".split('|')

  # Create a boolean mask where any cell in a row contains any of the substrings
  mask = df[columns_to_check].applymap(lambda x: contains_any_substring(x, substrings)).any(axis=1)

  interesting_rows = df[mask]

  mask = interesting_rows[columns_to_check].applymap(lambda x: contains_any_substring(x, unwanted_substrings)).any(axis=1)

  interesting_rows_reduced = interesting_rows[~mask]

  cols_to_keep = interesting_rows_reduced.columns
  columns_reduced = interesting_rows_reduced.loc[:, cols_to_keep]

  columns_reduced["MeaningfulDimension"] = interesting_rows_reduced["ConvertDimensions"].str.split('x', expand=True)[2].astype('int')

  file_size_describe = interesting_rows_reduced[["FileSize"]].describe()
  convertDimensions_describe = columns_reduced[["MeaningfulDimension"]].describe()

  wanted_size = 90 
  columns_reduced = columns_reduced.loc[columns_reduced['MeaningfulDimension'] >= wanted_size]

  columns_reduced["PathStart"] = columns_reduced["NiftiPath"].str.split('/',expand=True)[0]

  columns_reduced.to_csv('{}.csv'.format('uudet'), index = None, header=True)

  total_rows = columns_reduced.shape[0]
  print('total_rows',total_rows)

file_name = 'tags'
#get_interesting_rows(file_name)

