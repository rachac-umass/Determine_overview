# Folder and file names
data_folder = 'i2b2_data_files'

patient_file = 'i2b2_patient_details.parquet'
labresults_file = 'i2b2_labresults.parquet'
diag_file = 'i2b2_diagnoses.parquet'
med_file = 'i2b2_medications.parquet'
bmi_and_bp_file = 'i2b2_bmi_and_bp.parquet'
cvs_file = 'i2b2_cvs.parquet'


# Target columns names for each dataset
target_medication_columns = []
target_lab_results_columns = []
target_diag_columns= []
target_bmi_bp_columns = []
target_cvs_columns = []

# Column names required for final model
target_features = []


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

