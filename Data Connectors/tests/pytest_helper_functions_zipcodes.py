import pytest
from helper_functions import clean_zipcode

@pytest.fixture
def zipcodes():
    return {
        '5_digit': '12345',
        '5_dash_4': '12345-6789',
        '4_dash': '1234-5678',
        '4_no_dash': '1234',
        'invalid': '123',
        'empty': '',
        'non_numeric': 'abcd-efgh'
    }

def test_clean_zipcode_five_digits(zipcodes):
    assert clean_zipcode(zipcodes['5_digit']) == '12345'

def test_clean_zipcode_five_digits_with_dash(zipcodes):
    assert clean_zipcode(zipcodes['5_dash_4']) == '12345'

def test_clean_zipcode_four_digits_with_dash(zipcodes):
    assert clean_zipcode(zipcodes['4_dash']) == '01234'

def test_clean_zipcode_four_digits_no_dash(zipcodes):
    assert clean_zipcode(zipcodes['4_no_dash']) == '01234'

def test_clean_zipcode_invalid(zipcodes):
    assert clean_zipcode(zipcodes['invalid']) is None

def test_clean_zipcode_empty(zipcodes):
    assert clean_zipcode(zipcodes['empty']) is None

def test_clean_zipcode_non_numeric(zipcodes):
    assert clean_zipcode(zipcodes['non_numeric']) is None