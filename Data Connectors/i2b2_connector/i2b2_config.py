"""
Author: Chandra Harsha Rachabathuni
Date: 2025-10-23
Description: Python config file required for running "data_connector_i2b2.py".
Version: 1.0

TO-DO:
1. Add lab results and respective units for required features
2. Add units for weight, height for unit validation
"""

# Folder and file names
data_folder = 'i2b2_data_files'

patient_file = 'i2b2_patient_details.parquet'
labresults_file = 'i2b2_labresults.parquet'
diag_file = 'i2b2_diagnoses.parquet'
med_file = 'i2b2_medications.parquet'
bmi_and_bp_file = 'i2b2_bmi_and_bp.parquet'
cvs_file = 'i2b2_cvs.parquet'


# Lab results(loincs + required units) includes bmi and bp as well


# phemap paths
icd10_phe = './icd_code_to_phecode/Phecode_map_v1_2_icd10cm_beta.csv'
icd9_phe_v1 = './icd_code_to_phecode/phecode_map_v1_2_icd9.csv'
icd9_phe_v2 = './icd_code_to_phecode/phemap (1).rda'
icd9_to_icd10_mapping = './icd_code_to_phecode/icd9to10dictionary.txt'

#  Columns names required in each dataset
cohort_columns = ['PATIENT_NUM',
                   'Outcome',
                   'Gender_CD',
                   'Race_CD',
                   'Sex_CD',
                   'Hispanic_CD']

lab_results_columns = ['PATIENT_NUM',
                       'LabLOINC',
                       'Result_Number',
                       'Result_Unit']

diag_columns= ['PATIENT_NUM',
               'CONCEPT_CD']

medication_columns = ['PATIENT_NUM',
                      'rxcui']

bmi_bp_columns = ['PATIENT_NUM',
                  'CONCEPT_CD',
                  'NVAL_NUM',
                  'UNITS_CD']

cvs_columns = ['PATIENT_NUM',
               'Indicator',
               'FACT']

# Column names required for final model
target_features = ['Age_group', 'azithromycin', 'levothyroxine', 'acyclovir',
       'ceftriaxone', 'phe_401.1', 'phe_271.3', 'phe_41.0', 'phe_278.11',
       #'phe_649.1', 
        'LOINC:2085-9', 'LOINC:2345-7', 'LOINC:74774-1',
       'LOINC:27353-2', 'LOINC:9318-7', 'LOINC:62238-1', 'mode_height', 'BMI',
       'median_diastolic_value', 'ACS_MedHHIncome', 'ACS_pctCollGrad',
       'Race_CD_02', 'Race_CD_05', 'Hispanic_CD_Y', 'Gender_CD_M',
       'Gender_CD_W'] + ['LOINC:4548-4', 'LOINC:17856-6']


###### Other variables ######
ignore_rxnorm_codes = [
    "16681", "173", "1534763", "1368001", "1368384", "1368402", 
    "2627044", "1373458", "1545149", "2404", "1488564", "1486436", 
    "1727500", "1551291", "1545653", "1598392", "2281864", "1664314", 
    "1992672", "1992684", "1992825", "60548", "4816", "25789", 
    "647235", "606253", "4821", "352381", "25793", "4815", 
    "285129", "1007184", "51428", "1670007", "1727493", "139825", 
    "274783", "1858994", "400008", "631657", "1605101", "1008501", 
    "86009", "816726", "253182", "1100699", "1243019", "475968", 
    "1440051", "6809", "607999", "802646", "614348", "1043562", 
    "729717", "30009", "274332", "33738", "139953", "73044", 
    "84108", "857974", "1991302", "593411", "2638675", "2621880", 
    "2601723", "10633", "10635", "72610", "596554"
]

