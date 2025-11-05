WITH LatestVisit AS (
    SELECT DISTINCT PATIENT_NUM, Index_Start_date, Age_at_encounter_date
    FROM S348.dbo.DETERMINE_CRITERIA_IDX_DATE_AND_AGE
),

PositiveDates AS (
    -- Find patient id and dates that have a where condition satisfied
    SELECT 
        o.PATIENT_NUM,
        o.Start_date
    FROM 
        AAOCHIN2023.S348.OBSERVATION_FACT o
    WHERE 
        o.PATIENT_NUM IN (SELECT DISTINCT PATIENT_NUM FROM S348.dbo.DETERMINE_CRITERIA_IDX_DATE_AND_AGE)
        AND (
            o.CONCEPT_CD IN ( -- ICD 9 codes for type 2 diabetes
    'ICD9CM:250.00',
    'ICD9CM:250.02',
    'ICD9CM:250.10',
    'ICD9CM:250.12',
    'ICD9CM:250.20',
    'ICD9CM:250.22',
    'ICD9CM:250.30',
    'ICD9CM:250.32',
    'ICD9CM:250.40',
    'ICD9CM:250.42',
    'ICD9CM:250.50',
    'ICD9CM:250.52',
    'ICD9CM:250.60',
    'ICD9CM:250.62',
    'ICD9CM:250.70',
    'ICD9CM:250.72',
    'ICD9CM:250.80',
    'ICD9CM:250.82',
    'ICD9CM:250.90',
    'ICD9CM:250.92'
)
        OR (o.CONCEPT_CD LIKE 'ICD10CM:E11%') -- ICD 10 codes for type 2 diabetes
            -- LOINC codes for type 2 diabetes
        OR (o.NVAL_NUM < 20 AND o.NVAL_NUM >= 6.5 AND o.CONCEPT_CD = 'LOINC:4548-4' and o.UNITS_CD NOT IN ('mg/dL','mL/min'))
		OR (o.NVAL_NUM < 20 AND o.NVAL_NUM >= 6.5 AND o.CONCEPT_CD = 'LOINC:17856-6' and o.UNITS_CD NOT IN ('mg/dL'))
		OR (o.NVAL_NUM < 20 AND o.NVAL_NUM >= 6.5 AND o.CONCEPT_CD = 'LOINC:4549-2' and o.UNITS_CD NOT IN ('mg/dL','mL/min'))
		--OR ((o.NVAL_NUM < 20 AND o.NVAL_NUM >= 6.5 AND o.CONCEPT_CD = 'LOINC:41995-2' and o.UNITS_CD NOT IN ('mg/dL','mL/min'))
        )
),

FirstPositiveDate AS (
    SELECT 
        PATIENT_NUM, 
        MIN(Start_date) AS FirstOutcomeDate
    FROM 
        PositiveDates
    GROUP BY 
        PATIENT_NUM
),

EligibleWithPositive AS (
    SELECT  
        lv.PATIENT_NUM,
        lv.Index_Start_date,
        fp.FirstOutcomeDate,
        lv.Age_at_encounter_date,
        
        CASE 
            WHEN fp.FirstOutcomeDate IS NOT NULL 
                 AND fp.FirstOutcomeDate BETWEEN lv.Index_Start_date AND DATEADD(YEAR, 5, lv.Index_Start_date) THEN 1
            ELSE 0
        END AS Outcome
    FROM
        LatestVisit lv
    LEFT JOIN
        FirstPositiveDate fp ON lv.PATIENT_NUM = fp.PATIENT_NUM
    WHERE
        fp.FirstOutcomeDate IS NULL
        OR fp.FirstOutcomeDate > lv.Index_Start_date -- Exclusion of records which have OutcomeDate on or before Index_start_date
),

PATIENT_EXCLUDE_OTHER_DIABETES_BEFORE_VISIT_INDEX AS (
SELECT DISTINCT o.PATIENT_NUM 
FROM AAOCHIN2023.S348.OBSERVATION_FACT o 
LEFT JOIN
LatestVisit lv
ON o.PATIENT_NUM = lv.PATIENT_NUM
where
	o.PATIENT_NUM IN (
        SELECT DISTINCT PATIENT_NUM
    FROM S348.dbo.DETERMINE_CRITERIA_IDX_DATE_AND_AGE
    )
    AND ((o.CONCEPT_CD LIKE 'ICD10CM:E08%' OR 
 o.CONCEPT_CD LIKE 'ICD10CM:E09%' OR 
 --o.CONCEPT_CD LIKE 'ICD10CM:E10%' OR 
 o.CONCEPT_CD LIKE 'ICD10CM:E11%' OR 
 o.CONCEPT_CD LIKE 'ICD10CM:E13%' OR
 o.CONCEPT_CD LIKE 'ICD9CM:249%')
AND o.START_DATE <= lv.Index_Start_date)
),

PATIENT_EXLCUDE_TYPE10_GESTATIONAL AS (

SELECT DISTINCT o.PATIENT_NUM 
FROM AAOCHIN2023.S348.OBSERVATION_FACT o 
where
	o.PATIENT_NUM IN (
        SELECT DISTINCT PATIENT_NUM
    FROM S348.dbo.DETERMINE_CRITERIA_IDX_DATE_AND_AGE
    )
    AND (o.CONCEPT_CD LIKE 'ICD10CM:E10%' OR 
	o.CONCEPT_CD IN ('ICD9CM:250.01',
    'ICD9CM:250.03',
    'ICD9CM:250.11',
    'ICD9CM:250.13',
    'ICD9CM:250.21',
    'ICD9CM:250.23',
    'ICD9CM:250.31',
    'ICD9CM:250.33',
    'ICD9CM:250.41',
    'ICD9CM:250.43',
    'ICD9CM:250.51',
    'ICD9CM:250.53',
    'ICD9CM:250.61',
    'ICD9CM:250.63',
    'ICD9CM:250.71',
    'ICD9CM:250.73',
    'ICD9CM:250.81',
    'ICD9CM:250.83',
    'ICD9CM:250.91',
    'ICD9CM:250.93')
	OR o.CONCEPT_CD LIKE 'ICD9CM:648.0%'
	OR o.CONCEPT_CD LIKE 'ICD9CM:648.8%'
	OR o.CONCEPT_CD LIKE 'ICD10CM:O24.0%')
),


FinalCohort AS (
    SELECT 
        ewp.*,
        pd.DEATH_DATE,
        pd.Zip_CD
    FROM
        EligibleWithPositive AS ewp
    LEFT JOIN 
        AAOCHIN2023.S348.PATIENT_DIMENSION pd ON ewp.PATIENT_NUM = pd.PATIENT_NUM
    WHERE

        -- Exclusion of records which have death date before first positive outcome
        (pd.DEATH_DATE IS NULL OR pd.DEATH_DATE > ewp.FirstOutcomeDate)
        
        -- Exclude patients from the two subqueries
        AND ewp.PATIENT_NUM NOT IN (
            SELECT PATIENT_NUM FROM PATIENT_EXCLUDE_OTHER_DIABETES_BEFORE_VISIT_INDEX
        )
        AND ewp.PATIENT_NUM NOT IN (
            SELECT PATIENT_NUM FROM PATIENT_EXLCUDE_TYPE10_GESTATIONAL
        )
)


SELECT * INTO S348.dbo.DETERMINE_OCHIN_COHORT FROM FinalCohort;
