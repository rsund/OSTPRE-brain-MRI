import os
import glob
import bisect 
import xml.etree.ElementTree as ET

# Kaivaa xml tiedostosta IQR arvosanan ja laskee sen perusteella arvosanan matlab-kaavan mukaan
# mark2rps = @(mark) min(100,max(0,105 - real(mark)*10)) + isnan(real(mark)).*real(mark);  
# Tallentaa arvosanat grades.csv tiedostoon

def determine_grade(score):
    grades=['F', 'E-', 'E', 'E+', 'D-', 'D', 'D+', 'C-', 'C', 'C+', 'B-', 'B', 'B+', 'A-', 'A', 'A+']
    breakpoints=[]
    grade_range = 50/len(grades)
    for i in range(len(grades)):
        breakpoints.append(50+i*grade_range)
    i = bisect.bisect(breakpoints, score)
    return grades[i]

def check_grades(path):
    files = glob.glob(path + '/**/*.xml', recursive=True)
    total_files = len(files)
    for i, file in enumerate(files):
        print(f"{i}/{total_files}")
        is_cat_xml = 'cat_' in file
        if(not is_cat_xml):
            continue
        try:
            tree = ET.parse(file)
            root = tree.getroot()

            IQR_rating = root.find('qualityratings').find('IQR').text
            IQR_formatted = round(105-(float(IQR_rating)*10),2)
            grade = determine_grade(IQR_formatted)

            file_name = file.split('anat')[1].strip('/').strip('.xml')

            with open("grades.csv", 'a') as f:
                f.write(f"{grade};{IQR_formatted};{file_name}\n")

        except Exception as err:
            print(file)
            print('err', err)

cwd = os.getcwd()
print('cwd', cwd)
check_grades(cwd)
