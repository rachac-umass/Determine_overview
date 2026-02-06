"""
Author: Chandra Harsha Rachabathuni
Date: 2025-10-23
Description: Python file with helper functions required for running "data_connector_i2b2.py".
Version: 1.0

TO-DO:
1. Doc strings for remaining functions
"""

### importing required libraries ###
import polars as pl
import numpy as np
import requests
import pickle
import us
from tqdm import tqdm
import datetime
from scipy.stats import linregress
import os

from omop_config import *

### Helper functions ###

def check_columns(df: pl.DataFrame, target_columns: list[str])-> bool:
    """
    Check if the columns of the provided DataFrame match the target columns.

    Args:
        df (pd.DataFrame): The DataFrame to validate.
        target_columns (list[str]): The list of target column names for the DataFrame.
    """
    # is_subset = set(target_columns).issubset(df.columns)
    # print(target_columns)
    # print(df.columns)
    # print(is_subset)

    assert set(target_columns).issubset(df.collect_schema().names()), f'Target columns are not same as dataframe columns.'


def get_active_ingredient(rxcui):
       # Define the base URL of the RxNorm API
       base_url = "https://rxnav.nlm.nih.gov/REST/rxcui/"
       
       # Create the complete URL for the request
       request_url = f"{base_url}{rxcui}/related.json?tty=IN"
       
       try:
           # Make the GET request to the API
           response = requests.get(request_url)
           response.raise_for_status()  # Raise an error if the request was unsuccessful
        #    print(response)
           # Parse the JSON response
           data = response.json()

           # Extract the active ingredient from the JSON response
           ingredients = data.get('relatedGroup', {}).get('conceptGroup', [])
           for group in ingredients:
               if group.get('tty') == 'IN':  # Here we're looking for ingredients
                   active_ingredients = [concept['name'] for concept in group.get('conceptProperties', [])]
                   if active_ingredients == []:
                    # print(data)
                    continue
                   return active_ingredients
       
       except requests.RequestException as e:
           print(f"An error occurred: {e}")

       try:
            response = requests.get(f'https://rxnav.nlm.nih.gov/REST/rxcui/{rxcui}/historystatus.json?caller=RxNav')
            data = response.json()
                        
            active_ingredient_name = data['rxcuiStatusHistory']['definitionalFeatures'
                                                                   ]['ingredientAndStrength'
                                                                    ][0]['activeIngredientName']
            return active_ingredient_name
                
       except:
            return []
       

def meds_rxcui_to_api(list_rxcui, verbose = False):
    """
    Converts a list of RxNorm concept unique identifiers (RxCUIs) to their corresponding active ingredient information.

    This function uses a preloaded dictionary containing mappings from RxCUIs to active ingredients. 
    If an RxCUI doesn't have active ingredient data in the dictionary, it attempts to retrieve this information
    using the `get_active_ingredient` function.

    Args:
        list_rxcui (list): A list of RxNorm identifiers in "rx_cui:XXXXXX" format.
        verbose (bool): If True, prints RxCUIs for which active ingredient information could not be found. Default is False.

    Returns:
        list: A list of dataframes containing active ingredient information for each RxCUI.
              If no active ingredient is associated with a code, a placeholder list with the RxCUI is added.
    """
    # assert dict_rxnorm_active_ing is not None, "Parameter 'dict_rxnorm_active_ing' cannot be None"

    active_ingredients_list_dataframe = []
    no_act_ing_codes = []
    dict_rxnorm_active_ing = {}
    
    for rx_code in tqdm(list_rxcui):

        rx_code =str(rx_code)
        if ':' in rx_code:
            rx_code = rx_code.split(':')[-1].strip()

        try:
            if dict_rxnorm_active_ing[rx_code] == []:
                # print(rx_code)
                no_act_ing_codes.append(rx_code)
                active_ingredients_list_dataframe.append(['rx_cui:'+str(rx_code)]) ### handling codes that dont have active ingredient associated with them
                continue
            active_ingredients_list_dataframe.append(dict_rxnorm_active_ing[rx_code])
        except:
            active_ingredients_list_dataframe.append(get_active_ingredient(rx_code))
            dict_rxnorm_active_ing[rx_code] = get_active_ingredient(rx_code)

    if verbose:
        print("Rx code for which active ingredient cannot be found: ")
        print(no_act_ing_codes)

    print(dict_rxnorm_active_ing)
    return active_ingredients_list_dataframe

# def normalize_active_ingredients(name):
#     if type(name) == str:
#         return ('_').join(part.strip() for part in name.split('/'))
#     else:
#         print("value passed and its type: ",name, " and ", type(name))
#         raise ValueError("name parameter is not string")

def normalize_Active_ingredients(batch: pl.Series) -> pl.Series:
    # The function will be applied to a Series batch, so use .apply on the batch
    return batch.apply(
        lambda name: ('_').join(part.strip() for part in name.split('/'))
        if isinstance(name, str)
        else (
            print("value passed and its type: ", name, " and ", type(name)) or
            (_ for _ in ()).throw(ValueError("name parameter is not string"))
        )
    )

def normalize_active_ingredients_expr(col: pl.Expr) -> pl.Expr:
    # Use .str.split('/').list.join('_') for the string manipulation
    return col.str.split('/').list.join('_')


def read_icd_mappings(file_path):
    icd_mapping = {}

    # Open the file and read it line by line
    with open(file_path, 'r') as file:
        for line in file:
            # Strip any white space and split by the delimiter '|'
            parts = line.strip().split('|')
            if len(parts) >= 2:
                icd9_code = parts[0]
                icd10_code = parts[1]

                # Store the mapping in the dictionary
                icd_mapping[icd9_code] = icd10_code

    return icd_mapping


def get_phecode_from_concept_cd(cd, phewas_mapping_dicts):
    if cd.startswith('ICD10CM'):
        
        if cd.startswith("ICD10CM:Z"):
            return cd.split('.')[0]
        
        icd_code = cd.split(':')[-1]
        try:
            phecode = str(phewas_mapping_dicts['icd10_phe_dict'][icd_code])
            if phecode  =='1010.0':
                return cd.split('.')[0]
            return 'phe_'+phecode
        except:
            
            return cd.split('.')[0]
        
    elif cd.startswith('ICD9CM'):
        icd_code = cd.split(':')[-1]
        try:
            phecode = str(phewas_mapping_dicts['icd9_phe_dict_v2'][icd_code])
            if phecode =='1010.0':
                        return cd.split('.')[0]
            return 'phe_'+phecode
        except:
            try:
                phecode  = str(phewas_mapping_dicts['icd9_phe_dict_v1'][icd_code])
                if phecode =='1010.0':
                        return cd.split('.')[0]
                return 'phe_'+phecode
            
            except:
                try: #convert icd9 code to icd10 and try phemap
                    phecode = phewas_mapping_dicts['icd10_phe_dict'][phewas_mapping_dicts['icd9_to_10_mapping_dict'][icd_code]]
                    if phecode =='1010.0':
                        return cd.split('.')[0]
                    return 'phe_'+phecode
                except:
                    
                    return cd.split('.')[0]
    else:
        return cd
    

def convert_wt_kg_to_lb(kg: float)-> float:
    return kg * 2.20462

def clean_prefix_data(value: str)-> str:
    return value.split(':')[-1]


def custom_height_aggregator(ht_list):
    ht_list = [ht for ht in ht_list if (ht>36 and ht< 96)]
    if ht_list == []:
        return None
    return max(set(ht_list), key=ht_list.count)

def clean_zipcode(zipcode: str)-> str:
    if '-' in zipcode:
        zipcode = zipcode.split('-')[0]

    if len(zipcode)==4:
            return '0'+zipcode
    elif len(zipcode) ==5:
            return zipcode
    else:
            return None 
 


### SDoH extractor ###
def get_acs_data(census_object, zipstate_object, cds_fields, zipcode, year, missing_value = None):

    if zipcode is None:
        return_dict = dict(zip(cds_fields,[missing_value]*len(cds_fields)))
        return_dict['ACS_poverty'] = missing_value
        return_dict['ACS_unemployment'] = missing_value
        return_dict['pctCollGrad'] = missing_value

        return return_dict 
    
    
    try:
        state_fip = us.states.mapping('abbr', 'fips')[zipstate_object[zipcode].state]

        acs_data = census_object.acs5.state_zipcode(
        fields=cds_fields, 
        state_fips=state_fip,
        zcta=zipcode,
        year=year
        )

        # print(acs_data)

        try:
            acs_data[0]['ACS_poverty'] = np.round(((acs_data[0]['C17002_002E'] + acs_data[0]['C17002_003E'])/acs_data[0]['B01003_001E']) * 100,2)
        except:
            acs_data[0]['ACS_poverty'] =missing_value

        try:
            acs_data[0]['ACS_unemployment'] = np.round((acs_data[0]['B23025_005E']/acs_data[0]['B23025_001E'])*100,2)
        except:
            acs_data[0]['ACS_unemployment'] = missing_value
        try:
            acs_data[0]['pctCollGrad'] = np.round(((acs_data[0]['B15003_022E'] + acs_data[0]['B15003_023E'] + acs_data[0]['B15003_024E'] + acs_data[0]['B15003_025E'])/acs_data[0]['B01003_001E'])*100,2)
        except:
            acs_data[0]['pctCollGrad'] = missing_value

        return_dict = acs_data[0]




    except:
        #Handle invalid zipcodes
        return_dict = dict(zip(cds_fields,[missing_value]*len(cds_fields)))
        return_dict['ACS_poverty'] = missing_value
        return_dict['ACS_unemployment'] = missing_value
        return_dict['pctCollGrad'] = missing_value

    
    return return_dict


############# demographic functions #############
def func_map_race(value: str):

    if value is None:
        return "UNK"
    elif value not in race_map_dict.keys():
        return "UNK"
    else:
        return race_map_dict[value]
    
def func_map_sex(value: str):
    if value is None:
        return "UNK"
    elif value not in gender_map:
        return "UNK"
    else:
        return value
    
def func_map_ethinicity(value: str):
    if value is None:
        return "UNK"
    elif value not in ethinicity_map_dict.keys():
        return "UNK"
    else:
        return
    
def calculate_slope(subdf: pl.DataFrame) -> pl.DataFrame:
    dates = subdf["START_DATE"].to_numpy()
    values = subdf["Result_Number"].to_numpy()
    ordinal_dates = [d.astype(datetime.datetime).toordinal() for d in dates]
    
    if len(dates) > 3 and len(set(ordinal_dates)) > 1:
        slope, *_ = linregress(ordinal_dates, values)
    else:
        slope = 0.0

    # (You might want PATIENT_NUM in the return for identification)
    return pl.DataFrame({
        "PATIENT_NUM": [subdf["PATIENT_NUM"][0]],
        "slope": [slope]
    })


def loinc_codes_measurement_and_unit_check(lab_results_df, verbose, loinc_featureset_name = None , mode = None, loinc_features = []):

    lab_results_units_grouped_df = lab_results_df.group_by(['LabLOINC','Result_Unit']).agg([pl.col('Result_Number'),
                                                                                       pl.col('Result_Unit').count().alias('Unit_Count')
                                                                                       ]).sort(['LabLOINC','Result_Unit']).collect()
    if verbose:
        print("Number of unique combinations of loinc code and their units",len(lab_results_units_grouped_df))

    if loinc_features == []:
        raise ValueError("loinc_features are empty which is not valid input for the parameter")
    
    Loincs_with_more_than_1_unit_count = 0
    loinc_unit_pairs = []

    if loinc_featureset_name =='domain_expert':
        loincs_focus = loinc_features ## add code for obtaining de loinc codes

    elif loinc_featureset_name =='boruta_features':
        loincs_focus = loinc_features

    elif loinc_featureset_name =='all_features':
        loincs_focus = np.unique(lab_results_df.select('LabLOINC').collect()['LabLOINC'].to_list()).tolist()

    newpath = r'./logs/' 
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    
    with open(newpath + f'output_loincs_units_check_function_{loinc_featureset_name}_output.txt', 'w') as f:

        for loinc in loincs_focus:
            focus_group_loinc = lab_results_units_grouped_df.filter(pl.col('LabLOINC') == loinc).sort('Unit_Count', descending=True)
            
            if len(focus_group_loinc) == 1:
                continue
                
            total_count_patients_associated_with = focus_group_loinc['Unit_Count'].sum()
            if total_count_patients_associated_with < 1000:
                continue
            loinc_unit_with_max_count = focus_group_loinc['Result_Unit'][0]
            loinc_values_with_max_count = focus_group_loinc['Result_Number'][0].to_list()
            
            other_units_fg = focus_group_loinc.filter(pl.col('Result_Unit')!=loinc_unit_with_max_count)['Result_Unit'].to_list()
            max_val_in_pri_unit = np.quantile(loinc_values_with_max_count, 0.75)
            min_val_in_pri_unit = min(loinc_values_with_max_count)
            #q25 = np.quantile(loinc_values_with_max_count, 0.25)
            #q75 = np.quantile(loinc_values_with_max_count, 0.75)
            median_val_in_pri_unit = np.median(loinc_values_with_max_count)


        
            print_flag =  False
            for unit in other_units_fg:

                other_unit_values = focus_group_loinc.filter(pl.col('Result_Unit') == unit)['Result_Number'][0].to_list()
        #         print(other_unit_values)
                max_val_in_ou = max(other_unit_values)
                min_val_in_ou = min(other_unit_values)
                median_val_in_ou = np.mean(other_unit_values)

                if min_val_in_pri_unit <= median_val_in_ou <= max_val_in_pri_unit:
                    continue
                else:
                    loinc_unit_pairs.append((loinc,unit))
                    f.write("Unit to focus: ",unit)            
                    Loincs_with_more_than_1_unit_count += 1
                    if median_val_in_pri_unit/median_val_in_ou > 10 or median_val_in_pri_unit/median_val_in_ou < 0.1:
                        f.write(" The value is orders above 10th order\n")
                        print_flag = True
                    else:
                        f.write(" The value is orders below 10th order\n")
                        print_flag = True

            if print_flag:
                f.write(str(focus_group_loinc))

            f.write('*'*20 +"\n\n")

        f.write(f"\n\n\nNumber of unique loinc codes with disparities {len(loinc_unit_pairs)}")

        f.close()
    
    if mode == 'standardize_remove_outlier':
        return loinc_unit_pairs

    
        

