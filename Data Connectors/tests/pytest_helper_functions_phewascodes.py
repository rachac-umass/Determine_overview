import pytest

from helper_functions import meds_rxcui_to_api




######################## Active  ingredient ########################
# Mock function for get_active_ingredient
def mock_get_active_ingredient(rx_code):
    mock_data = {
        "123456": [{"ingredient": "Active1"}],
        "654321": [{"ingredient": "Active2"}],
    }
    return mock_data.get(rx_code, [])

@pytest.fixture
def rxcui_dict():
    return {
        "123456": [{"ingredient": "Active1"}],  # valid entry
        "234567": [],  # no active ingredient entry
        "345678": [{"ingredient": "Active3"}],  # valid entry
    }

def test_meds_rxcui_to_api_successful_conversion(rxcui_dict, monkeypatch):
    monkeypatch.setattr(__builtins__, 'get_active_ingredient', mock_get_active_ingredient)
    
    list_rxcui = ["rxcui:123456", "rxcui:345678"]
    
    result = meds_rxcui_to_api(list_rxcui, rxcui_dict)
    
    assert result == [
        [{"ingredient": "Active1"}],
        [{"ingredient": "Active3"}]
    ], "Test failed for successful conversion of RxCUIs."

def test_meds_rxcui_to_api_no_active_ingredient(rxcui_dict, monkeypatch):
    monkeypatch.setattr(__builtins__, 'get_active_ingredient', mock_get_active_ingredient)
    
    list_rxcui = ["rxcui:234567"]
    
    result = meds_rxcui_to_api(list_rxcui, rxcui_dict)
    
    assert result == [
        ["rxcui:234567"]
    ], "Test failed for RxCUIs without active ingredients."

def test_meds_rxcui_to_api_fetch_from_api(rxcui_dict, monkeypatch):
    monkeypatch.setattr(__builtins__, 'get_active_ingredient', mock_get_active_ingredient)
    
    list_rxcui = ["rxcui:654321"]
    
    result = meds_rxcui_to_api(list_rxcui, rxcui_dict)
    
    assert result == [
        [{"ingredient": "Active2"}]
    ], "Test failed for RxCUIs requiring API fetch."

def test_meds_rxcui_to_api_empty_list(rxcui_dict):
    list_rxcui = []
    
    result = meds_rxcui_to_api(list_rxcui, rxcui_dict)
    
    assert result == [], "Test failed for empty RxCUI list."

def test_meds_rxcui_to_api_verbose_flag(rxcui_dict, monkeypatch, capsys):
    monkeypatch.setattr(__builtins__, 'get_active_ingredient', mock_get_active_ingredient)
    
    list_rxcui = ["rxcui:234567", "rxcui:789012"]
    
    meds_rxcui_to_api(list_rxcui, rxcui_dict, verbose=True)
    
    captured = capsys.readouterr()
    assert "Rx code for which active ingrident cannot be found:" in captured.out


######################## Phewas code ########################


