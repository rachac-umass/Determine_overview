"""
Author: Chandra Harsha Rachabathuni
Date: 2025-10-23
Description: Python script to transform omop raw data into data required for type-2 diabetes prediction model.
Version: 1.0

TO-DO:
1. Add data validator function in helper functions
2. Verbose functionality
3. Skip feature engineering
4. Move varibles to omop config (there are lots of variables to move)
5. Scaling features according to model type to be used
"""

import omop_config
from helper_functions import check_columns,  meds_rxcui_to_api, get_active_ingredient
from helper_functions import read_icd_mappings, get_initial_char_for_icd_codes
from helper_functions import clean_prefix_data
from helper_functions import custom_height_aggregator

import argparse
import polars as pl
import numpy as np
import pyreadr


### Command line aruments ###
parse = argparse.ArgumentParser(description = '''This script is for transforming raw omop data into dataset that can be used for model prediction
                                
                                Data requirements:
                                    Folder that contains all required files:
                                        1. patient info
                                        2. medications
                                        3. lab results
                                        4. diagnoses
                                        5. bmi and bp
                                        6. community vital signs
                                
                                *** Make sure the filenames are modified to match the names present in the config file. ***
                                *** Column names used for model prediction are in config file. ***
                                *** Handles both csv and parquet files. Recommended: parquet ***

                                ''',
                                epilog = 'Refer to Github repo Wiki for data requirements')

parse.add_arguments('--data_folder',
                    type = str, 
                    help = 'Path to the folder that contains all required input files.')

parse.add_arguments('--use_lab_results_aggregation', 
                    action = argparse.BooleanOptionalAction, 
                    default = True, 
                    help = 'Argument that allows lab results aggregation(units filter and median agg). By default, it is True')

parse.add_arguments('--medication_tranformation',
                    type = str,
                    choices = ['rxcui_api','api'], ### add atc later 
                    help = '''Argument to to choose whether to transform the medications into type.
                    Options:
                    1. rxcui -> rxcui to active ingredient
                    2. active ingredient -> no transformation
                    3. ATC -> rxcui to ATC''')


parse.add_arguments('--file_format',
                   type = str,
                   default = 'parquet',
                   choices = ['csv','parquet'],
                   help = 'Argument to choose file format of input data format. Defualt value is parquet(pl)')

parse.add_arguments('--use_bmi_bp_aggregation',
                   action = argparse.BooleanOptionalAction, 
                    default = True, 
                    help = 'Argument that allows bmi and bp aggregation(median agg). By default, it is True')


parse.add_arguments('--output_before_joining',
                   action = argparse.BooleanOptionalAction, 
                    default = True, 
                    help = 'Argument that allows to write transformed datasets before joining to local')

parse.add_arguments('--output_file_format',
                   type = str,
                   default = 'parquet',
                    choices = ['parquet','cvs'], 
                    help = 'Argument that allows to change output file format(csv or parquet)')


parse.add_arguments('-v','--verbose',action = argparse.BooleanOptionalAction, default = False,
                    help = 'Enable verbose mode (optional)')

args = parse.parse_args()

class script_config:
    # folder and files to use #
    data_folder = args.data_folder

    cohort_file = data_folder + omop_config.patient_file
    lab_file = data_folder + omop_config.labresults_file
    diag_file = data_folder + omop_config.diag_file
    med_file = data_folder + omop_config.med_file
    bmi_bp_file = data_folder + omop_config.bmi_and_bp_file
    cvs_file = data_folder + omop_config.cvs_file

    use_lab_results_aggregation = False

    # verbose
    verbose = args.ver

    # data format of input files
    file_format = args.file_format

    # column name checks #
    patient_file_columns = omop_config.patient_columns

    # Medications transformation argument 
    medications_transformation = args.medication_tranformation

    # Final data columns required #
    final_columns_required = omop_config.modeling_columns


############# Main code to run #################
if __name__ == '__main__':

    ####### Loading all require files(Features and data for helper functions) #######
    if script_config.file_format == 'pl':
        cohort_df = pl.scan_parquet(script_config.patient_file)
        lab_df = pl.scan_parquet(script_config.labresults_file)
        diag_df = pl.scan_parquet(script_config.diag_file)
        meds_df = pl.scan_parquet(script_config.med_file)
        bmi_bp_df = pl.scan_parquet(script_config.bmi_and_bp_file)
        cvs_df = pl.scan_parquet(script_config.cvs_file)
    else:
        cohort_df = pl.scan_csv(script_config.patient_file)
        lab_df = pl.scan_csv(script_config.labresults_file)
        diag_df = pl.scan_csv(script_config.diag_file)
        meds_df = pl.scan_csv(script_config.med_file)
        bmi_bp_df = pl.scan_csv(script_config.bmi_and_bp_file)
        cvs_df = pl.scan_csv(script_config.cvs_file)

    # Diagnoses mapping functions
    phewas_mapping_dicts = {"icd10_phe" : pl.read_csv('./icd_code_to_phecode/Phecode_map_v1_2_icd10cm_beta.csv'),
               "icd9_phe_v1" :  pl.read_csv('./icd_code_to_phecode/phecode_map_v1_2_icd9.csv',infer_schema_length=100000),
               "icd9_phe_v2" : pyreadr.read_r("./icd_code_to_phecode/phemap (1).rda")['phemap'].drop_duplicates(subset=['icd9']),
                "icd9_to_10_mapping_dict" : read_icd_mappings('./icd_code_to_phecode/icd9to10dictionary.txt') 
    }

    # Required lab results for modeling
    lab_loinc_and_unit_tuple = []


    #### Column names validation ####
    check_columns(cohort_df, omop_config.cohort_columns)
    check_columns(lab_df, omop_config.lab_results_columns)
    check_columns(diag_df, omop_config.diag_columns)
    check_columns(meds_df, omop_config.medication_columns)
    check_columns(bmi_bp_df, omop_config.bmi_bp_columns)
    check_columns(cvs_df, omop_config.cvs_columns)

    ##### Required columns from each dataset for modeling ####
    target_medication_columns = omop_config.target_medication_columns
    target_lab_results_columns = omop_config.target_lab_results_columns
    target_diag_columns = omop_config.target_diag_columns
    target_bmi_bp_columns = omop_config.target_bmi_bp_columns
    target_cvs_columns = omop_config.target_cvs_columns

    ############################### Data transformation for required columns ############################### 

    ############## Transform medications ##############
    if script_config.medications_transformation == 'api':
        assert 'Active Ingredient' in meds_df.columns
        pass

    elif script_config.medications_transformation == 'rxcui_api':
        meds_df = meds_df.with_columns(pl.col['rxcui'].map_batches(meds_rxcui_to_api).alias('Active_ingrident'))

    # Filtering patients who are on diabetes type 2 medications
    act_list_to_drop = [get_active_ingredient(rxcui) for rxcui in omop_config.ignore_rxnorm_codes]
    medications_to_drop = ["_".join(sublist) if len(sublist) > 1 else sublist[0] for sublist in act_list_to_drop]
    med_ignore_patient_num = np.unique(meds_df.filter(pl.col('Active_ingrident').is_in(medications_to_drop)).collect()['PATIENT_NUM'].to_list())

    ############## Transform diagnoses codes to PheWas codes/ATC codes ##############
    diag_df = diag_df.with_columns(pl.col("phecode_map").map_elements(get_initial_char_for_icd_codes, return_dtype = pl.Utf8).alias("phecode_map"))

    ############## Transform lab results: Filter required lab results and unit standardization and aggregation (median) ##############

    # Define keys for the dictionaries
    keys = ("LabLOINC", "Unit")

    # Convert tuple to list of dictionaries
    dict_loinc_unit = [dict(zip(keys, item)) for item in lab_loinc_and_unit_tuple]

    lab_df = lab_df.filter(
        pl.struct(["LabLOINC", "Units"]).is_in(dict_loinc_unit)
    )
    lab_df = lab_df.group_by(['PATIENT_NUM','LabLOINC']).agg(pl.col('Result_Number').median().alias('Result_Number'))


    ############## Transform bmi and bp ##############

    # *** BMI *** 
    wt_df = bmi_bp_df.filter(pl.col('CONCEPT_CD').str.contains('WT'))
    wt_df = wt_df.with_columns(pl.col("NVAL_NUM").cast(pl.Float32))
    wt_df = wt_df.filter((pl.col('NVAL_NUM') > 65) & (pl.col('NVAL_NUM') < 600))
    average_weight_per_patient = wt_df.group_by("PATIENT_NUM").agg([
        pl.col("NVAL_NUM").median().alias("median_weight")
        ])
    
    height_df = bmi_bp_df.filter(pl.col('CONCEPT_CD').str.contains('WT'))
    
    mode_height_df = height_df.group_by("PATIENT_NUM").agg([
    pl.col("NVAL_NUM").map_elements(custom_height_aggregator, return_dtype = pl.Float32).alias("mode_height")
    ])

    # Join the DataFrames on 'patient_id'
    bmi_ht_wt_df = mode_height_df.join(average_weight_per_patient, on="PATIENT_NUM", how="inner")

    # Calculate BMI
    bmi_ht_wt_df = bmi_ht_wt_df.with_columns(
        (pl.col("median_weight") / (pl.col("mode_height") ** 2) * 703).alias("BMI")
    )

    # ** Diastolic and systolic BP **
    dia_bp_df = bmi_bp_df.filter(pl.col('CONCEPT_CD').str.contains('VIT\|DIA'))
    sys_bp_df = bmi_bp_df.filter(pl.col('CONCEPT_CD').str.contains('VIT\|SYS'))
    dia_bp_df = bmi_bp_df.with_columns(pl.col('NVAL_NUM').cast(pl.Float32).alias('NVAL_NUM'))
    sys_bp_df = bmi_bp_df.with_columns(pl.col('NVAL_NUM').cast(pl.Float32).alias('NVAL_NUM'))

    dia_bp_df  = dia_bp_df.filter((pl.col('NVAL_NUM') > 40) &
                             (pl.col('NVAL_NUM') < 160))
    sys_bp_df  = sys_bp_df.filter((pl.col('NVAL_NUM') > 60) &
                                (pl.col('NVAL_NUM') < 190))
    
    average_dia_bp_per_patient = dia_bp_df.group_by("PATIENT_NUM").agg([
        pl.col("NVAL_NUM").mean().alias("median_diastolic_value"),

    ])

    average_sys_bp_per_patient = sys_bp_df.group_by("PATIENT_NUM").agg([
        pl.col("NVAL_NUM").mean().alias("median_systolic_value"),

    ])

    ############## Transform patient/cohort file ##############

    cohort_df = cohort_df.with_columns(pl.col('Age_at_encounter_date').map_elements(lambda x : '18-34' if x <= 34 \
                                                        else '35-44' if 35<=x<=44 else '45-54'\
                                                        if 45<=x<=54 else '55-64' if 55<=x<=64 else \
                                                        '65-74' if  65<=x<=74 else '75_older').alias('Age_group'))
    
    cohort_df = cohort_df.with_columns( pl.col('Sex_CD').map_elements(clean_prefix_data, return_dtype = pl.Utf8).alias('Sex_CD'),
                                               pl.col('Race_CD').map_elements(clean_prefix_data, return_dtype = pl.Utf8).alias('Race_CD'),
                                               pl.col('Hispanic_CD').map_elements(clean_prefix_data, return_dtype = pl.Utf8).alias('Hispanic_CD'),
                                               pl.col('Gender_CD').map_elements(clean_prefix_data, return_dtype = pl.Utf8).alias('Gender_CD')
                                  )
    
    ## Race  ##
    cohort_df = cohort_df.with_columns(pl.col('Race_CD'))

    # Values to be replaced
    values_to_replace = ["06", "07"]
    replacement_value = "UN"

    cohort_df = cohort_df.with_columns([
    pl.col('Sex_CD').str.replace(r'^(NI|OT|UN)$', 'UNK', literal=False).alias('Sex_CD'),
    pl.col('Race_CD').str.replace(r'^(NI|OT|UN)$', 'UNK', literal=False).alias('Race_CD'),
    pl.col('Hispanic_CD').str.replace(r'^(R|NI|UN)$', 'UNK', literal=False).alias('Hispanic_CD'),
    pl.col('Gender_CD').str.replace(r'^(OT|NI|UN)$', 'UNK', literal=False).alias('Gender_CD')
    ])

    # Replace multiple values with a single value
    cohort_df = cohort_df.with_columns(
        pl.when(pl.col("Race_CD").is_in(values_to_replace))
        .then(pl.lit(replacement_value))
        .otherwise(pl.col("Race_CD"))
        .alias("Race_CD")
    )

    age_group_dict = {'18-34': 0,
                  '35-44': 1,
                  '45-54':2,
                  '55-64':3,
                  '65-74':4,
                  '75_older':5    
    }
    cohort_df = cohort_df.with_columns(pl.col('Age_group').replace_strict(age_group_dict).alias('Age_group'))


    ####### Pivoting and joining the datasets #######
    meds_df = meds_df.with_columns(pl.lit(1).alias('usage'))
    medications_pivot_df = meds_df.pivot(on = 'Active_ingrident', index = 'PATIENT_NUM', values = 'usage', aggregate_function = 'max')

    lab_results_pivot_df = lab_df.pivot(on = 'LabLOINC', index = 'PATIENT_NUM', values = 'Result_Number')
    
    diag_df = diag_df.with_columns(pl.lit(1).alias('Usage'))
    diagnoses_pivot_df = diag_df.pivot(on = 'phecode_map', index = 'PATIENT_NUM', values = 'Usage', aggregate_function = 'max')

    community_vital_pivot_df = cvs_df.pivot(on = 'Indicator', 
                                                    index ='PATIENT_NUM', 
                                                    values ='FACT', 
                                                    aggregate_function = 'mean')
    community_vital_pivot_df = community_vital_pivot_df.fill_null(-100)


    ###### Removing patients with active ingredient of diabetes from required patient ids #####
    med_pids = np.unique(meds_df.select(['PATIENT_NUM']).collect()['PATIENT_NUM'].to_list())
    dx_pids = np.unique(diag_df.select(['PATIENT_NUM']).collect()['PATIENT_NUM'].to_list())
    lab_pid = np.unique(lab_df.select(['PATIENT_NUM']).collect()['PATIENT_NUM'].to_list())
    req_patients = set(med_pids) | set(dx_pids) | set(lab_pid)

    req_patients -= set(med_ignore_patient_num)

    ##################### Joining the pivoted feature dataframes #####################
    cohort_df = cohort_df.filter(pl.col('PATIENT_NUM').is_in(req_patients))

    temp_df = cohort_df.join(medications_pivot_df, on = 'PATIENT_NUM', how = 'left')
    temp_df = temp_df.fill_null(0)

    temp_df = temp_df.join(diagnoses_pivot_df, on = 'PATIENT_NUM', how = 'left')
    temp_df = temp_df.fill_null(0)

    temp_df = temp_df.join(lab_results_pivot_df, on = 'PATIENT_NUM', how = 'left')
    temp_df = temp_df.fill_null(-100)

    

    temp_df = temp_df.join(bmi_ht_wt_df, on = 'PATIENT_NUM', how = 'left')
    temp_df = temp_df.fill_null(-100)
    temp_df = temp_df.join(dia_bp_df, on = 'PATIENT_NUM', how = 'left')
    temp_df = temp_df.fill_null(-100)
    temp_df = temp_df.join(sys_bp_df, on = 'PATIENT_NUM', how = 'left')
    temp_df = temp_df.fill_null(-100)

    temp_df = temp_df.join(cvs_df, on = 'PATIENT_NUM',how ='left')
    temp_df = temp_df.fill_null(-100)

    # Define the columns to one-hot encode
    columns_to_one_hot_encode = ['Sex_CD', 'Race_CD', 'Hispanic_CD', 'Gender_CD']

    # Combine the non-encoded columns with the encoded ones
    columns_to_one_hot_encode = ['Sex_CD', 'Race_CD', 'Hispanic_CD', 'Gender_CD']

    df_final_modeing_dataset = temp_df.to_dummies(columns=columns_to_one_hot_encode)

    df_final_modeing_dataset = df_final_modeing_dataset.drop('Sex_CD_UNK','Race_CD_UNK','Hispanic_CD_UNK','Gender_CD_UNK')

    ########################## Filtering only features required for modeling ###########################
    df_final_modeing_dataset = df_final_modeing_dataset.filter(pl.col(omop_config.target_features))

