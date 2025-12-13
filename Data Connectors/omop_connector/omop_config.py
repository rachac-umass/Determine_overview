"""
Author: Chandra Harsha Rachabathuni
Date: 2025-10-23
Description: Python config file required for running "data_connector_omop.py".
Version: 1.0

TO-DO:
"""

# Folder and file names
data_folder = 'omop_data_files'

patient_file = 'omop_patient_details'
labresults_file = 'omop_labresults'
diag_file = 'omop_diagnoses'
med_file = 'omop_medications'
bmi_and_bp_file = 'omop_bmi_and_bp'
cvs_file = 'omop_cvs'

# secret key and fields for census data
census_key = ''
cvs_fields = ["B19083_001E", # GINI INDEX OF INCOME INEQUALITY
            "B19013_001E", # Medican household income
            "C17002_002E", # poverty <50
            "C17002_003E", # poverty 50<x<99
            "B17001_002E", # Total population
            "B15002_001E", # SEX BY EDUCATIONAL ATTAINMENT FOR THE POPULATION 25 YEARS AND OVER
             "B01003_001E",
          "B23025_001E", # EMPLOYMENT STATUS FOR THE POPULATION 16 YEARS AND OVER (TOTAL)
          "B23025_005E", #EMPLOYMENT STATUS FOR THE POPULATION 16 YEARS AND OVER (UNEMPLOYED)
          "B15003_022E",
          "B15003_023E",
          "B15003_024E",
          "B15003_025E"
         ]

missing_numerical_replace_val = -10

# Lab results(loincs + required units) includes bmi and bp as well

# phemap paths
icd10_phe = './icd_code_to_phecode/Phecode_map_v1_2_icd10cm_beta.csv'
icd9_phe_v1 = './icd_code_to_phecode/phecode_map_v1_2_icd9.csv'
icd9_phe_v2 = './icd_code_to_phecode/phemap (1).rda'
icd9_to_icd10_mapping = './icd_code_to_phecode/icd9to10dictionary.txt'


#  Columns names required in each dataset
patient_columns = ['person_id',
                  'Age_at_index'
                   'Outcome',
                   'Race_CD',
                   'Hispanic_CD',
                   'Sex_CD',
                #    'Encounter',
                   'Zipcode'] ### Need Gender, Race and Hispanic columns

lab_results_columns = ['person_id',
                       'measurement_source_value', # LOINC code
                       'value_source_value', # lab result measurement value
                       'unit_source_value'] # Unit of measurement

diag_columns= ['person_id',
               'vocabulary_id', ## ICD9 or ICD10
               'condition_source_value'] ## ICD code or concept_code

medication_columns = ['person_id',
                      'concept_code', ## rxnorm code from table concept
                      'concept_name', ## Active ingrident
                      ]

cvs_columns = ['person_id',
               'ACS_GINI',
                'ACS MedHHIncome',
                'ACS pctPoverty100',
                'ACS_Unemployment',
                'ACS_pctCollGrad'] #needs extraction from aqhs

# Column names required for final model or modeling dataset
target_features = ['Age_group', 'azithromycin', 'levothyroxine', 'acyclovir',
       'ceftriaxone', 'phe_401.1', 'phe_271.3', 'phe_41.0', 'phe_278.11',
       #'phe_649.1', 
        'LOINC:2085-9', 'LOINC:2345-7', 'LOINC:74774-1',
       'LOINC:27353-2', 'LOINC:9318-7', 'LOINC:62238-1', 'mode_height', 'BMI',
       'median_diastolic_value', 'ACS_MedHHIncome', 'ACS_pctCollGrad',
       'Race_CD_02', 'Race_CD_05', 'Hispanic_CD_Y', 'Gender_CD_M',
       'Gender_CD_W'] + ['LOINC:4548-4', 'LOINC:17856-6']

target_label = ['Outcome']

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

weight_loincs = [
    "29463-7",
    "3141-9",
    "3142-7",
    "75292-3",
    "79348-9",
    "8350-1",
    "8351-9"
]

weight_loinc_unit = ['kg','lb']

height_loincs = [
    "3137-7",
    "3138-5",
    "8302-2",
    "8306-3",
    "8308-9"
]

height_loinc_unit = ['in']

systolic_loinc_codes = [
    "75997-7",
    "76215-3",
    "76534-7",
    "8479-8",
    "8480-6",
    "8508-4"
]

diastolic_loinc_codes = [
    "75995-1",
    "76213-8",
    "76535-4",
    "8462-4",
    "8496-2"
]

dia_sys_stolic_unit = ['mm Hg','mm[Hg]', 'mmHg']

diastolic_range = [40,160]
systolic_range = [60,190]

bmi_range = [15,60]
height_range = [36,96]
weight_range = [65,600]

#### Demogrpahics dict ####
race_map_dict = {
    "American Indian or Alaska Native": "01",
    "Asian": "02",
    "Black or African American": "03",
    "Native Hawaiian or Other Pacific Islander": "04",
    "White": "05",
}

ethinicity_map_dict  = {
    "Hispanic or Latino": "Y",
    "Not Hispanic or Latino": "N",
}

gender_map = ['F','M']

acs_Cvs_fields = fields = ["B19083_001E", # GINI INDEX OF INCOME INEQUALITY
            "B19013_001E", # Medican household income
            "C17002_002E", # poverty <50
            "C17002_003E", # poverty 50<x<99
            "B17001_002E", # Total population 
            "B15002_001E", # SEX BY EDUCATIONAL ATTAINMENT FOR THE POPULATION 25 YEARS AND OVER 
             "B01003_001E",
          "B23025_001E", # EMPLOYMENT STATUS FOR THE POPULATION 16 YEARS AND OVER (TOTAL)
          "B23025_005E", #EMPLOYMENT STATUS FOR THE POPULATION 16 YEARS AND OVER (UNEMPLOYED)
          "B15003_022E",
          "B15003_023E",
          "B15003_024E",
          "B15003_025E"
         ]
census_key = '' #### KEY GOES HERE
############################################## UNITS VALIDATION DICT ##############################################

units_validation_tuple_boruta = {
    "LOINC:2085-9": "mg/dL", #[100 -1000]
    "LOINC:2345-7": "mg/dL",
    "LOINC:74774-1": "mg/dL",
    "LOINC:27353-2": "mg/dL",
    "LOINC:9318-7": "ug/mg{creat}",
    "LOINC:62238-1": "mL/min/{1.73_m2}"
}