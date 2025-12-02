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

from omop_config import *

### Helper functions ###

def check_columns(df: pl.DataFrame, target_columns: list[str])-> bool:
    """
    Check if the columns of the provided DataFrame match the target columns.

    Args:
        df (pd.DataFrame): The DataFrame to validate.
        target_columns (list[str]): The list of target column names for the DataFrame.
    """
    assert target_columns in df.columns , f'Target columns are not same as dataframe columns.'


def get_active_ingredient(rxcui):
       # Define the base URL of the RxNorm API
       base_url = "https://rxnav.nlm.nih.gov/REST/rxcui/"
       
       # Create the complete URL for the request
       request_url = f"{base_url}{rxcui}/related.json?tty=IN"
       
       try:
           # Make the GET request to the API
           response = requests.get(request_url)
           response.raise_for_status()  # Raise an error if the request was unsuccessful
           print(response)
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
            response = requests.url(f'https://rxnav.nlm.nih.gov/REST/rxcui/{rxcui}/historystatus.json?caller=RxNav')
            data = response.json()
                        
            active_ingredient_name = data['rxcuiStatusHistory']['definitionalFeatures'
                                                                   ]['ingredientAndStrength'
                                                                    ][0]['activeIngredientName']
            return active_ingredient_name
                
       except:
            return []
       

def meds_rxcui_to_api(list_rxcui, dict_rxnorm_active_ing = None, verbose = False):
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
    assert dict_rxnorm_active_ing is not None, "Parameter 'dict_rxnorm_active_ing' cannot be None"

    active_ingredients_list_dataframe = []
    no_act_ing_codes = []
    
    for rxnorm in list_rxcui:
        rx_code = rxnorm.split(':')[-1].strip()
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
        print("Rx code for which active ingrident cannot be found: ")
        print(no_act_ing_codes)

    return active_ingredients_list_dataframe

def normalize_active_ingridents(name):
    return ('_').join(name.split(','))


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
    

def convert_wt_kg_to_lb(kg):
    return kg * 2.20462

def clean_prefix_data(value):
    return value.split(':')[-1]


def custom_height_aggregator(ht_list):
    ht_list = [ht for ht in ht_list if (ht>36 and ht< 96)]
    if ht_list == []:
        return None
    return max(set(ht_list), key=ht_list.count)



### SDoH extractor ###
def get_acs_data(census_object, zipstate_object, cds_fields, zipcode, year, missing_value = None):

    invalid_zipcodes = []
    #Handling zipcode with "-". Example: 01610-3301
    if '-' in zipcode:
        zipcode = zipcode.split('-')[0]
          
    try:
        state_fip = us.states.mapping('abbr', 'fips')[zipstate_object[zipcode].state]
    
        acs_data = census_object.acs5.state_zipcode(
        fields=cds_fields, 
        state_fips=state_fip,
        zcta=zipcode,
        year=year
        )

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
            acs_data[0]['ACS_unemployment'] = missing_value

        return_dict = acs_data[0]
    except:
        # Handle invalid zipcodes
        return_dict = dict(zip(cds_fields,[missing_value]*len(cds_fields)))
    
    return return_dict


############# demogrpahics function #############
def func_map_race(value):

    if value is None:
        return "UNK"
    elif value not in race_map_dict.keys():
        return "UNK"
    else:
        return race_map_dict[value]
    
def func_map_sex(value):
    if value is None:
        return "UNK"
    elif value not in gender_map_dict.keys():
        return "UNK"
    else:
        return 
    
def func_map_ethinicity(value):
    if value is None:
        return "UNK"
    elif value not in ethinicity_map_dict.keys():
        return "UNK"
    else:
        return