import pytest
import numpy as np

# import the function to test
# from your_module import get_acs_data
from helper_functions import get_acs_data
import us
from pyzipcode import ZipCodeDatabase
from census import Census

@pytest.fixture
def cds_fields():
    return [
        "B19083_001E", # GINI INDEX OF INCOME INEQUALITY
        "B19013_001E", # Median household income
        "C17002_002E",
        "C17002_003E", # EMPLOYMENT STATUS FOR THE POPULATION 16 YEARS AND OVER
        "B17001_002E",
        "B01003_001E",
    ]

@pytest.fixture
def zipcode():
    return '12345'

@pytest.fixture
def zipstate_object(zipcode):
    # Just needs a .state attribute
    class ZipState:
        state = 'CA'
    return {zipcode: ZipState()}

@pytest.fixture
def state_fip():
    return '06'

@pytest.fixture
def census_object():
    class DummyACS5:
        def state_zipcode(self, fields, state_fips, zcta, year):
            # Simulate dictionary output
            return [{
                # arbitrary values for each needed census field
                "B19083_001E": 20,
                "B19013_001E": 40000,
                "C17002_002E": 100,
                "C17002_003E": 50,
                "B17001_002E": 150,
                "B01003_001E": 1000,
                # more for unemployment
                "B23025_005E": 50, # Unemp numerator
                "B23025_001E": 1000, # Unemp denom
                # college grad fields
                "B15003_022E": 10,
                "B15003_023E": 15,
                "B15003_024E": 25,
                "B15003_025E": 50,
            }]
    class DummyCensus:
        acs5 = DummyACS5()
    return DummyCensus()
    
@pytest.fixture(autouse=True)
def patch_us_states(monkeypatch, state_fip):
    class DummyUsStates:
        @staticmethod
        def mapping(src, dst):
            # abbr to fips
            return {'CA': state_fip}
    import sys
    sys.modules['us'] = __import__('types') # ensure 'us' import works
    setattr(sys.modules['us'], 'states', DummyUsStates())

def test_get_acs_data_zipcode_none(census_object, zipstate_object, cds_fields):

    result = get_acs_data(
        census_object,
        zipstate_object,
        cds_fields,
        zipcode=None,
        year=2020,
        missing_value=-1
    )
    # All values should be set to missing_value in result
    for k in cds_fields + ['ACS_poverty', 'ACS_unemployment', 'pctCollGrad']:
        assert result[k] == -1

def test_get_acs_data_happy_path(census_object, zipstate_object, cds_fields, zipcode):

    result = get_acs_data(
        census_object,
        zipstate_object,
        cds_fields,
        zipcode,
        year=2021,
        missing_value=np.nan
    )
    # Field values
    assert result["B19083_001E"] == 20
    assert result["B19013_001E"] == 40000
    # Poverty: ((100+50)/1000) * 100 = 15.0
    assert np.isclose(result['ACS_poverty'], 15.0)
    # Unemployment: (50/1000)*100 = 5.0
    assert np.isclose(result['ACS_unemployment'], 5.0)
    # College: ((10+15+25+50)/1000)*100 = 10.0
    assert np.isclose(result['pctCollGrad'], 10.0)

def test_get_acs_data_exception(monkeypatch, census_object, zipstate_object, cds_fields, zipcode):
    # Patch census_object.acs5.state_zipcode to raise
    class DummyACS5Fail:
        def state_zipcode(self, **kwargs):
            raise Exception("fail")
    census_object.acs5 = DummyACS5Fail()

    result = get_acs_data(
        census_object,
        zipstate_object,
        cds_fields,
        zipcode,
        year=2021,
        missing_value="fail"
    )
    for k in cds_fields + ['ACS_poverty', 'ACS_unemployment', 'pctCollGrad']:
        assert result[k] == "fail"