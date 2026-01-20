"""
Author: Chandra Harsha Rachabathuni
Date: 2025-10-23
Description: Python script to transform omop raw data into data required for type-2 diabetes prediction model.
Version: 1.0

TO-DO:
1. Add data validator function in helper functions (units validator for sys and dia bp)
2. Verbose functionality
3. Skip feature engineering

5. Scaling features according to model type to be used
"""

import omop_config
from helper_functions import check_columns,  meds_rxcui_to_api, get_active_ingredient, normalize_active_ingridents, normalize_active_ingredients_expr
from helper_functions import read_icd_mappings, get_phecode_from_concept_cd
from helper_functions import custom_height_aggregator,convert_wt_kg_to_lb
from helper_functions import func_map_race, func_map_sex, func_map_ethinicity
from helper_functions import clean_zipcode, get_acs_data

import argparse
import polars as pl
import numpy as np
import pyreadr
from pyzipcode import ZipCodeDatabase

from census import Census
from us import states

import sys


### Command line aruments ###
parse = argparse.ArgumentParser(
    description="""
This script transforms raw OMOP data into a dataset suitable for model prediction.

Data Requirements:
  - Folder containing all required files:
      1. Patient info
      2. Medications
      3. Lab results
      4. Diagnoses
      5. BMI and BP
      6. Community vital signs

*** Ensure filenames match those specified in the config file. ***
*** Model prediction columns are defined in the config file. ***
*** Handles both CSV and Parquet files (Parquet recommended). ***
    """,
    epilog='Refer to the GitHub repo Wiki for detailed data requirements.',
    formatter_class=argparse.RawTextHelpFormatter
)

parse.add_argument(
    '--data_folder',
    type=str,
    required=True,
    help='Path to the input folder containing all required OMOP files (e.g., "data/input/")'
)

parse.add_argument(
    '--use_lab_results_aggregation', 
    action=argparse.BooleanOptionalAction, 
    default=True,
    help='Aggregate lab results by median after unit harmonization. Recommended: keep enabled. Default: %(default)s'
)

parse.add_argument(
    '--medication_transformation',
    type=str,
    choices=['rxcui_api', 'api'],
    help=(
        "Method for medication transformation:\n"
        "  - 'rxcui_api': converts RxCUI to active ingredient (default).\n"
        "  - 'api': treats RxCUI as active ingredient (no transformation).\n"
        "  (ATC will be supported later.)"
    )
)

parse.add_argument(
    '--input_file_format',
    type=str,
    default='parquet',
    choices=['csv', 'parquet'],
    help='File format for input data. Choose between "csv" and "parquet". Default: %(default)s'
)

parse.add_argument(
    '--use_bmi_bp_aggregation',
    action=argparse.BooleanOptionalAction,
    default=True,
    help='Aggregate BMI and blood pressure by median. Default: %(default)s'
)

parse.add_argument(
    '--output_before_joining',
    action=argparse.BooleanOptionalAction,
    default=False,
    help=(
        "If enabled, saves individual transformed datasets to disk "
        "before merging them for modeling or prediction. Default: %(default)s"
    )
)

parse.add_argument(
    '--missing_numerical_value_negative_10',
    action=argparse.BooleanOptionalAction,
    default=False,
    help=(
        "If enabled, sets missing values for numerical columns to -10. "
        "Recommended: keep enabled. Default: %(default)s"
    )
)

parse.add_argument(
    '--retrieve_sdoh_cvs',
    action=argparse.BooleanOptionalAction,
    default=False,
    help=(
        "If enabled, collects SDoH zipcode-level from ACS."
        "Use this if cvs_data file is not available. Default: %(default)s"
    )
)

parse.add_argument(
    '--apply_unit_standardizer_lab_results',
    action=argparse.BooleanOptionalAction,
    default=False,
    help=(
        "If enabled, does unit standardization by removing invalid units."
        "Default: %(default)s, meaning script uses lab results from boruta features and it's corresponding units."
    )
)


parse.add_argument(
    '--output_file_format',
    type=str,
    default='parquet',
    choices=['parquet', 'csv'],
    help='Output file format for the processed dataset. Choose "csv" or "parquet". Default: %(default)s'
)

parse.add_argument(
    '--mode',
    type=str,
    default='predict',
    choices=['predict', 'model'],
    help='Choose the operation mode: "model" (for training), "predict" (for inference). Default: %(default)s'
)

parse.add_argument(
    '--feature_set',
    type=str,
    default='Boruta',
    choices=['Boruta', 'All'],
    help=(
        'Features to include in the final processed dataset:\n'
        '  - "Boruta": Only features selected by the Boruta algorithm.\n'
        '  - "All": All available features (no selection). Default: %(default)s'
    )
)

parse.add_argument(
    '--verbose',
    action=argparse.BooleanOptionalAction,
    default=False,
    help=(
        "If enabled, outputs progress in the command line at each step of the script. "
        "Default: %(default)s."
    )
)

parse.add_argument(
    '--census_api_key',
    type = str,
    default = None,
    help=(
        "(Required) Your U.S. Census API key. This argument must be provided to access Census data."
    )
)


args = parse.parse_args()

class script_config:
    # folder and files to use #
    data_folder = args.data_folder

    patient_file = data_folder + omop_config.patient_file
    labresults_file = data_folder + omop_config.labresults_file
    diag_file = data_folder + omop_config.diag_file
    med_file = data_folder + omop_config.med_file
    bmi_bp_file = data_folder + omop_config.bmi_and_bp_file
    cvs_file = data_folder + omop_config.cvs_file

    use_lab_results_aggregation = False

    # verbose
    verbose = args.verbose

    # data format of input files
    file_format = args.input_file_format

    # column name checks #
    patient_file_columns = omop_config.patient_columns

    # Medications transformation argument 
    medication_transformation = args.medication_transformation

    # Final data columns required #
    final_columns_required = omop_config.target_modeling_features


############# Main code to run #################
if __name__ == '__main__':

    ####### Loading all require files(Features and data for helper functions) #######
    if script_config.verbose:
        print("Loading the required dataset files.")
    if script_config.file_format == 'pl':
        suffix = '.pl'
        cohort_df = pl.scan_parquet(script_config.patient_file + suffix)
        lab_df = pl.scan_parquet(script_config.labresults_file+ suffix)
        diag_df = pl.scan_parquet(script_config.diag_file+ suffix)
        meds_df = pl.scan_parquet(script_config.med_file+ suffix)
        if args.retrieve_sdoh_cvs:
            cvs_df = pl.scan_parquet(script_config.cvs_file+ suffix)
    else:
        suffix ='.csv'
        cohort_df = pl.scan_csv(script_config.patient_file+ suffix)
        lab_df = pl.scan_csv(script_config.labresults_file+ suffix)
        diag_df = pl.scan_csv(script_config.diag_file+ suffix)
        meds_df = pl.scan_csv(script_config.med_file+ suffix)
        if not args.retrieve_sdoh_cvs:
            cvs_df = pl.scan_csv(script_config.cvs_file+ suffix)


    if script_config.verbose:
        print("Loading the files required for mapping icd codes to phewas.")

    # Diagnoses mapping functions
    phewas_mapping_dicts = {"icd10_phe" : pl.read_csv(omop_config.icd10_phe),
               "icd9_phe_v1" :  pl.read_csv(omop_config.icd9_phe_v1, infer_schema_length=100000),
               "icd9_phe_v2" : pyreadr.read_r(omop_config.icd9_phe_v2)['phemap'].drop_duplicates(subset=['icd9']),
                "icd9_to_10_mapping_dict" : read_icd_mappings(omop_config.icd9_to_icd10_mapping) 
    }


    if script_config.verbose:
        print("Column validation across all files(cohort, lab results, diagnoses, medications and cvs)")
    #### Column names validation ####
    check_columns(cohort_df, omop_config.patient_columns)
    check_columns(lab_df, omop_config.lab_results_columns)
    check_columns(diag_df, omop_config.diag_columns)
    check_columns(meds_df, omop_config.medication_columns)
    # check_columns(bmi_bp_df, omop_config.bmi_bp_columns) #### Not required for omop as bmi and bp are extracted from lab results
    if not args.retrieve_sdoh_cvs:
        check_columns(cvs_df, omop_config.cvs_columns)

    ##### Required columns from each dataset for modeling ####
    # target_medication_columns = omop_config.target_medication_columns
    # target_lab_results_columns = omop_config.target_lab_results_columns
    # target_diag_columns = omop_config.target_diag_columns
    # target_bmi_bp_columns = omop_config.target_bmi_bp_columns
    # if not args.retrieve_sdoh_cvs:
    #     target_cvs_columns = omop_config.target_cvs_columns

    ############################### Data transformation for required columns ############################### 

    ############## Transform medications ##############
    if script_config.medication_transformation == 'api':
        assert 'ancestor_drug_concept_name' in meds_df.collect_schema().names()
        # meds_df = meds_df.with_columns(pl.col('ancestor_drug_concept_name').map_batches(normalize_active_ingridents, return_dtype = pl.Utf8).alias('Active_ingrident'))
        meds_df = meds_df.with_columns(
                        normalize_active_ingredients_expr(pl.col('ancestor_drug_concept_name')).alias('Active_ingrident')
                        )

    elif script_config.medication_transformation == 'rxcui_api':
        rxcuis_list = meds_df.select(pl.col('concept_code')).collect()['concept_code'].to_list()
        act_ing_list = meds_rxcui_to_api(rxcuis_list)

        final_act_ing_list = ["_".join(sublist) if len(sublist) > 1 else sublist[0] for sublist in act_ing_list]

        meds_df = meds_df.with_columns(Active_ingrident = pl.Series(final_act_ing_list))

    # Filtering patients who are on diabetes type 2 medications
    act_list_to_drop = [get_active_ingredient(rxcui) for rxcui in omop_config.ignore_rxnorm_codes]
    medications_to_drop = ["_".join(sublist) if len(sublist) > 1 else sublist[0] for sublist in act_list_to_drop]
    med_ignore_patient_id = np.unique(meds_df.filter(pl.col('Active_ingrident').is_in(medications_to_drop)).collect()['person_id'].to_list())

    ############## Transform diagnoses codes to PheWas codes/ATC codes ##############
    diag_df = diag_df.with_columns((pl.col("vocabulary_id")+':'+pl.col('concept_code')).alias('ICD_code'))
    diag_df = diag_df.with_columns(pl.col("ICD_code").map_elements(get_phecode_from_concept_cd, return_dtype = pl.Utf8).alias("phecode_map"))

    ############## Transform lab results: Filter required lab results and unit standardization and aggregation (median) ##############

    # Define keys for the dictionaries
    keys = ("measurement_source_value", "unit_source_value")

    # renaming columns
    lab_df = lab_df.rename({"measurement_source_value": "LabLOINC", "unit_source_value": "Units","value_as_number":'Result_Number'})

    # Convert tuple to list of dictionaries
    # dict_loinc_unit = [dict(zip(keys, item)) for item in omop_config.lab_loinc_and_unit_tuple]

    #### Filtering boruta features ####
    # lab_final_df = lab_df.filter(
    #     pl.struct(["LabLOINC", "Units"]).is_in(dict_loinc_unit)
    # )

    combo_df = pl.DataFrame(tuple(omop_config.units_validation_tuple_boruta.items()), schema=['LabLOINC', 'Result_Unit'])
    
    # Create a struct column in combo_df for unique combinations
    combo_structs = combo_df.select(pl.struct(["LabLOINC", "Result_Unit"]).alias("combo_struct"))["combo_struct"]

    # Now filter lab_results_df to exclude rows with LabLOINC/Result_Unit combos in combo_df
    lab_final_df = lab_df.filter(
        pl.struct(["LabLOINC", "Result_Unit"]).is_in(combo_structs)
    )

    lab_final_df = lab_final_df.group_by(['person_id','LabLOINC']).agg(pl.col('Result_Number').median().alias('Result_Number'))


    ############## Transform bmi and bp ##############

    # *** BMI *** 
    wt_df = lab_df.filter(pl.col('LabLOINC').is_in(omop_config.weight_loincs))
    assert wt_df.select('Units').collect()['Units'].is_in(omop_config.weight_loinc_unit).all(), "weight df units column contains values in other than 'kg' or 'lb'"
    
    ### Uncomment below later
    
    # wt_df =wt_df.with_columns(
        # pl.when(pl.col('Units')=='kg').then(pl.col('Result_Number').map_elements(convert_wt_kg_to_lb)).otherwise(pl.col('Result_Number')).alias('Result_Number')
    # )
    wt_df = wt_df.with_columns(pl.col("Result_Number").cast(pl.Float32))
    wt_df = wt_df.filter((pl.col('Result_Number') > 65) & (pl.col('Result_Number') < 600))
    average_weight_per_patient = wt_df.group_by("person_id").agg([
        pl.col("Result_Number").median().alias("median_weight")
        ])
    
    height_df = lab_df.filter(pl.col('LabLOINC').is_in(omop_config.height_loincs))
    
    # Asserting units be in inches
    assert height_df.select('Units').collect()['Units'].is_in(omop_config.height_loinc_unit).all(), "height df units column contains values other than 'inches' or 'in'"
    
    mode_height_df = height_df.group_by("person_id").agg([
    pl.col("Result_Number").map_elements(custom_height_aggregator, return_dtype = pl.Float32).alias("mode_height")
    ])

    # Join the DataFrames on 'person_id'
    bmi_ht_wt_df = mode_height_df.join(average_weight_per_patient, on="person_id", how="inner")

    # Calculate BMI
    bmi_ht_wt_df = bmi_ht_wt_df.with_columns(
        (pl.col("median_weight") / (pl.col("mode_height") ** 2) * 703).alias("BMI")
    )

    # print("bmi_wht_wt_df")
    # print(bmi_ht_wt_df.collect().head())

    bmi_ht_wt_df = bmi_ht_wt_df.filter((pl.col('BMI') > omop_config.bmi_range[0]) &
                             (pl.col('BMI') < omop_config.bmi_range[1]))

    # ** Diastolic and systolic BP **
    dia_bp_df = lab_df.filter(pl.col('LabLOINC').is_in(omop_config.diastolic_loinc_codes))
    sys_bp_df = lab_df.filter(pl.col('LabLOINC').is_in(omop_config.systolic_loinc_codes))
    dia_bp_df = dia_bp_df.with_columns(pl.col('Result_Number').cast(pl.Float32).alias('Result_Number'))
    sys_bp_df = sys_bp_df.with_columns(pl.col('Result_Number').cast(pl.Float32).alias('Result_Number'))


    #### Add units validation for dia and sys bp
    assert dia_bp_df.select('Units').collect()['Units'].is_in(omop_config.dia_sys_stolic_unit).all(), "diastolic df units column contains values other than 'mm Hg', 'mmHg' or 'mm[Hg]'"
    assert sys_bp_df.select('Units').collect()['Units'].is_in(omop_config.dia_sys_stolic_unit).all(), "systolic df units column contains values other than 'mm Hg', 'mmHg' or 'mm[Hg]'"


    dia_bp_df  = dia_bp_df.filter((pl.col('Result_Number') > omop_config.diastolic_range[0]) &
                             (pl.col('Result_Number') < omop_config.diastolic_range[1]))
    sys_bp_df  = sys_bp_df.filter((pl.col('Result_Number') > omop_config.systolic_range[0]) &
                                (pl.col('Result_Number') < omop_config.systolic_range[1]))
    
    average_dia_bp_per_patient = dia_bp_df.group_by("person_id").agg([
        pl.col("Result_Number").mean().alias("median_diastolic_value"),

    ])

    average_sys_bp_per_patient = sys_bp_df.group_by("person_id").agg([
        pl.col("Result_Number").mean().alias("median_systolic_value"),

    ])

    # print(average_dia_bp_per_patient.fetch(5))
    # print(average_dia_bp_per_patient.fetch(5))
    # print(bmi_ht_wt_df.collect().head())

    ############## Transform patient/cohort file (demographics) ##############
    cohort_df = cohort_df.with_columns(pl.col('Age_at_index').map_elements(lambda x : '18-34' if x <= 34 \
                                                        else '35-44' if 35<=x<=44 else '45-54'\
                                                        if 45<=x<=54 else '55-64' if 55<=x<=64 else \
                                                        '65-74' if  65<=x<=74 else '75_older').alias('Age_group'))
    


    
    cohort_df = cohort_df.with_columns( pl.col('Sex_CD').map_elements(func_map_sex, return_dtype = pl.Utf8).alias('Sex_CD'),
                                               pl.col('Race_CD').map_elements(func_map_race, return_dtype = pl.Utf8).alias('Race_CD'),
                                               pl.col('Hispanic_CD').map_elements(func_map_ethinicity, return_dtype = pl.Utf8).alias('Hispanic_CD'),
                                            #    pl.col('Gender_CD').map_elements(clean_prefix_data, return_dtype = pl.Utf8).alias('Gender_CD')
                                  )


    age_group_dict = {'18-34': 0,
                  '35-44': 1,
                  '45-54':2,
                  '55-64':3,
                  '65-74':4,
                  '75_older':5    
    }
    cohort_df = cohort_df.with_columns(pl.col('Age_group').replace_strict(age_group_dict).alias('Age_group'))
    print(cohort_df.collect()['Race_CD'])


    ############## Creating SDoH CVS ##############
    if args.retrieve_sdoh_cvs:

        zcdb = ZipCodeDatabase()
        census = Census(args.census_api_key)

        try:
            # Minimal call: Get the name and total population (B01003_001E) for Washington state (FIPS 53)
            # from the 2020 Decennial Census SF1 dataset
            test_data = census.acs5.get(('NAME', 'B25034_010E'), {'for': 'state:{}'.format(states.MD.fips)})

            if test_data and isinstance(test_data, list) and len(test_data) > 0:
                print("API key is active and working.")
                print(f"Successfully retrieved test data: {test_data}")
            else:
                print("API key may be inactive or invalid. Received no data.")

        except Exception as e:
            print(f"An error occurred: {e}")
            print("API key is likely invalid or inactive.")
            print("Ensure your key is activated via the email link provided by the Census Bureau.")
            sys.exit(1)

        

        pat_zipcodes_data = cohort_df.select('Zipcode').collect()['Zipcode'].to_list()

        pat_zipcodes_data_clean = [clean_zipcode(str(zip)) for zip in pat_zipcodes_data]
        pat_zipcodes_data_clean_uni = np.unique(pat_zipcodes_data_clean)
        acs_data_zips = {}
                    
        for pat_zip in pat_zipcodes_data_clean_uni:

            acs_data_zips[pat_zip] = get_acs_data(census, zcdb, 
                                          omop_config.cvs_fields,
                                          pat_zip, 
                                          2017)## Add code to dynamically add zipcode
            
        pat_acs_data = []
        for zip in pat_zipcodes_data_clean:
            pat_acs_data.append(acs_data_zips[zip])

        cvs_df = pl.DataFrame(pat_acs_data).lazy()

        cvs_df = cvs_df.rename({'B19083_001E': 'ACS_GINI',
                                  'B19013_001E': 'ACS_MedHHIncome',
                                   'ACS_poverty':'ACS_pctPoverty100',
                                    'ACS_unemployment': 'ACS_Unemployment',
                                    'pctCollGrad': 'ACS_pctCollGrad'})
        

    ####### Pivoting and joining the datasets #######
    meds_df = meds_df.with_columns(pl.lit(1).alias('usage'))
    medications_pivot_df = meds_df.pivot(on = 'Active_ingrident', index = 'person_id', values = 'usage', aggregate_function = 'max')

    lab_results_pivot_df = lab_df.pivot(on = 'LabLOINC', index = 'person_id', values = 'Result_Number')
    
    diag_df = diag_df.with_columns(pl.lit(1).alias('Usage'))
    diagnoses_pivot_df = diag_df.pivot(on = 'phecode_map', index = 'person_id', values = 'Usage', aggregate_function = 'max')

    community_vital_pivot_df = cvs_df#.pivot(on = 'Indicator', 
                                                    # index ='person_id', 
                                                    # values ='FACT', 
                                                    # aggregate_function = 'first')
    
    #######
    
    if args.missing_numerical_value_negative_10:
        community_vital_pivot_df = community_vital_pivot_df.fill_null(omop_config.missing_numerical_replace_val)
    # community_vital_pivot_df = community_vital_pivot_df.fill_null(-100)


    ###### Removing patients with active ingredient of diabetes from required patient ids #####
    med_pids = np.unique(meds_df.select(['person_id']).collect()['person_id'].to_list())
    dx_pids = np.unique(diag_df.select(['person_id']).collect()['person_id'].to_list())
    lab_pid = np.unique(lab_df.select(['person_id']).collect()['person_id'].to_list())
    req_patients = set(med_pids) | set(dx_pids) | set(lab_pid)

    req_patients -= set(med_ignore_patient_id)

    ##################### Joining the pivoted feature dataframes #####################
    cohort_df = cohort_df.filter(pl.col('person_id').is_in(req_patients))

    temp_df = cohort_df.join(medications_pivot_df, on = 'person_id', how = 'left')
    temp_df = temp_df.fill_null(0)

    temp_df = temp_df.join(diagnoses_pivot_df, on = 'person_id', how = 'left')
    temp_df = temp_df.fill_null(0)

    temp_df = temp_df.join(lab_results_pivot_df, on = 'person_id', how = 'left')
    if args.missing_numerical_value_negative_10:
        temp_df = temp_df.fill_null(omop_config.missing_numerical_replace_val)

    temp_df = temp_df.join(bmi_ht_wt_df, on = 'person_id', how = 'left')
    if args.missing_numerical_value_negative_10:
        temp_df = temp_df.fill_null(omop_config.missing_numerical_replace_val)
    temp_df = temp_df.join(dia_bp_df, on = 'person_id', how = 'left')
    if args.missing_numerical_value_negative_10:
        temp_df = temp_df.fill_null(omop_config.missing_numerical_replace_val)
    temp_df = temp_df.join(sys_bp_df, on = 'person_id', how = 'left')
    if args.missing_numerical_value_negative_10:
        temp_df = temp_df.fill_null(omop_config.missing_numerical_replace_val)

    temp_df = temp_df.join(cvs_df, on = 'person_id',how ='left')
    if args.missing_numerical_value_negative_10:
        temp_df = temp_df.fill_null(omop_config.missing_numerical_replace_val)

    # Define the columns to one-hot encode
    columns_to_one_hot_encode = ['Sex_CD', 'Race_CD', 'Hispanic_CD']#, 'Gender_CD']

    df_final_modeing_dataset = temp_df.to_dummies(columns=columns_to_one_hot_encode)

    df_final_modeing_dataset = df_final_modeing_dataset.drop('Sex_CD_UNK','Race_CD_UNK','Hispanic_CD_UNK')#,'Gender_CD_UNK')

    ########################## Filtering only features required for modeling ###########################
    if args.mode == 'predict':
        if args.feature_set == 'Boruta':
            df_final_modeing_dataset = df_final_modeing_dataset.filter(pl.col(omop_config.target_features))
    elif args.mode == 'evaluate':
        if args.feature_set == 'Boruta':
            df_final_modeing_dataset = df_final_modeing_dataset.filter(pl.col(omop_config.target_features + omop_config.target_label))

