from bids_validator import BIDSValidator
import pandas as pd
import shutil
import os
import csv

# Luo BIDS-kansiorakenteen ja tarkista onko se oikein

def validate(path):
    return BIDSValidator().is_bids(path)

def create_folders(file_name):
  df = pd.DataFrame(pd.read_csv('{}.csv'.format(file_name)))
  
  # TODO
  # Luo headerin, tämän voisi tehdä järkevämmin
  with open('participants.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['participant_id', 'sex', 'age'])
  
  for i, row in enumerate(df.iterrows()):
      print(i)
      path = row[1]['NiftiPath']
      sub = row[1]['PathStart']
      patient_sex = row[1]['PatientSex']
      
      year = int(row[1]['PatientBirthDate'].split('-')[0])
      month = int(row[1]['PatientBirthDate'].split('-')[1])
      day = int(row[1]['PatientBirthDate'].split('-')[2])

      acc_year = int(row[1]['AcquisitionDateTime'].split('-')[0])
      acc_month = int(row[1]['AcquisitionDateTime'].split('-')[1])
      acc_day = int(row[1]['AcquisitionDateTime'].split('-')[2].split('T')[0]) #Esimerkki aikaleimasta: "2012-05-14T10:48:16.390000"

      # TODO
      # Tämän voisi tehdä hieman selkeämmin
      month_over = (month - acc_month) > 0
      day_over = (day - acc_day) > 0

      has_birthday_after_image = 0
      if(month_over and day_over):
        has_birthday_after_image = -1

      patient_age = acc_year-year + has_birthday_after_image

      session = path.split('/')[1].replace('0', '')
      file_name = path.split('/')[2]

      start = session
      end = '.nii.gz'
      old_file_name = file_name[file_name.find(start)+len(start):file_name.rfind(end)]
      try: 
          top_folder = 'bids-folders'
          path_to_create = os.path.join(top_folder, sub, session, 'anat') 

          os.makedirs(path_to_create)   

          # T1w -kuvia tallennetaan
          new_path = "{}/{}_{}".format(path_to_create, sub, file_name.replace(old_file_name, '_T1w'))
          is_valid = validate(new_path)
          if(not is_valid):
             print(new_path, 'is not valid')
             exit(99)

          new_path_json = "{}/{}_{}".format(path_to_create, sub, file_name.replace('nii.gz', 'json').replace(old_file_name,'_T1w'))
          new_path_log = "{}/{}_{}".format(path_to_create, sub, 'log.txt')

          json_path = 'nifti/' + path.replace('nii.gz', 'json')
          formatted_path = 'nifti/' + path

          shutil.copyfile(formatted_path, new_path)
          shutil.copyfile(json_path, new_path_json)
          f = open(new_path_log, "a")
          f.write(old_file_name)
          f.close()

          # Kirjoitetaan tiedot participants.tsv tiedostoon
          # alkuperäinen tiedostopolku, sukupuoli ja ikä
          with open('participants.tsv', 'a') as out_file:
              tsv_writer = csv.writer(out_file, delimiter='\t')
              tsv_writer.writerow([sub, patient_sex, patient_age])

      except OSError as error:
          print(error)  
  df.to_csv('{}.csv'.format('validated'), index = None, header=True)
  print('Done')


file_name = 'columns_reduced'
create_folders(file_name)