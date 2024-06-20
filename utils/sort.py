import pandas as pd

# Muuttaa xlsx tiedoston csv:ksi
# Kerää kuvista kiinnostavat, ja tallentaa ne csv:ksi
# csv:n perusteella voidaan luoda BIDS-kansiorakenne

def convert_xlsx_to_csv(file_name):
  excel_data_df = pd.read_excel('{}.xlsx'.format(file_name))
  excel_data_df.to_csv('{}.csv'.format(file_name), 
                    index = None,
                    header=True)

def get_interesting_rows(file_name):
  df = pd.DataFrame(pd.read_csv('{}.csv'.format(file_name)))

  # Halutut protokollanimet: "mpr", "mprage", "t1_se", 'sag', 'tra', 'T1', 'TFE'
  # Poistettavat protokollanimet: "gd", "GD", "Gd", "gD"
  wanted_words = "mpr|mprage|t1_se|sag|tra|T1|TFE"
  unwanted_words = "gd|GD|Gd|gD"
  column_to_check = "ProtocolName" 

  interesting_rows = df[df[column_to_check].str.contains(wanted_words)==True]

  interesting_rows_reduced = interesting_rows[~interesting_rows[column_to_check].str.contains(unwanted_words)==True]

  cols_to_keep = ['NiftiPath', 'FileSize', 'ConvertDimensions', 'PatientBirthDate', 'PatientSex', 'AcquisitionDateTime']
  columns_reduced = interesting_rows_reduced.loc[:, cols_to_keep]

  # Dimensioiden perusteella voidaan päätellä onko kuva kokonainen
  columns_reduced["MeaningfulDimension"] = columns_reduced["ConvertDimensions"].str.split('x', expand=True)[2].astype('int')

  # 256x256x155x1
  file_size_describe = interesting_rows_reduced[["FileSize"]].describe()
  convertDimensions_describe = columns_reduced[["MeaningfulDimension"]].describe()

  print('Possibly interesting info:')
  print(file_size_describe)
  print(convertDimensions_describe)

  wanted_size = 90 
  columns_reduced = columns_reduced.loc[columns_reduced['MeaningfulDimension'] >= wanted_size]

  columns_reduced["PathStart"] = columns_reduced["NiftiPath"].str.split('/',expand=True)[0]
  columns_reduced = columns_reduced.groupby(['PathStart','FileSize','NiftiPath', 'MeaningfulDimension', 'PatientBirthDate','PatientSex', 'AcquisitionDateTime']).size().reset_index(name='count')

  columns_reduced.to_csv('{}.csv'.format('columns_reduced'), index = None, header=True)

  deliveryTypes = df["ManufacturersModelName"].value_counts()
  print('Valmistajat ja määrät: \n', deliveryTypes)

# Konvertoidaan excel csv:ksi, koska csv on helpompi käsitellä ja vie vähemmän tilaa
# convert_xlsx_to_csv('tags')

file_name = 'tags'
get_interesting_rows(file_name)

