import pandas as pd
import pyodbc
import sqlalchemy as db
from sqlalchemy import text
from sqlalchemy import create_engine, inspect
from sqlalchemy.engine import URL

import argparse
import os

def get_alchemy_engine(dbtype, server, database, uid=None, pwd=None, port=None, driver=None):
    dbtype_lower = dbtype.lower()

    # ── SQL Server / MSSQL ────────────────────────────────────────────────────
    if dbtype_lower in ('mssql', 'sqlserver', 'sql_server'):
        drivers = pyodbc.drivers()
        if not driver:
            for recommended in [
                "ODBC Driver 18 for SQL Server",
                "ODBC Driver 17 for SQL Server",
                "SQL Server",
            ]:
                if recommended in drivers:
                    driver = recommended
                    break
            else:
                raise ValueError('No suitable SQL Server ODBC driver found')
        else:
            if driver not in drivers:
                raise ValueError(
                    f"Specified driver '{driver}' not in available drivers: {drivers}"
                )
        print(f"Using driver: {driver}")

        if uid and pwd:
            con_str = (
                fr'Driver={driver};'
                fr'Server={server};'
                fr'Database={database};'
                fr'UID={uid};'
                fr'PWD={pwd};'
                r'TrustServerCertificate=Yes;'
            )
        else:
            # Windows (Kerberos/SSPI) auth
            con_str = (
                fr'Driver={driver};'
                fr'Server={server};'
                fr'Database={database};'
                r'Trusted_Connection=yes;'
                r'TrustServerCertificate=Yes;'
            )

        conn_url = URL.create("mssql+pyodbc", query={"odbc_connect": con_str})
        engine = create_engine(conn_url)
        return engine

    # ── MariaDB ───────────────────────────────────────────────────────────────
    elif dbtype_lower in ('mariadb', 'maria_db'):
        if not (uid and pwd):
            raise ValueError("MariaDB requires uid and pwd; Windows auth is not supported.")
        port = port or 3306
        conn_url = URL.create(
            "mariadb+mariadbconnector",
            username=uid,
            password=pwd,
            host=server,
            port=port,
            database=database,
        )
        engine = create_engine(conn_url)
        return engine

    # ── MySQL ─────────────────────────────────────────────────────────────────
    elif dbtype_lower in ('mysql', 'my_sql'):
        if not (uid and pwd):
            raise ValueError("MySQL requires uid and pwd.")
        port = port or 3306
        conn_url = URL.create(
            "mysql+pymysql",
            username=uid,
            password=pwd,
            host=server,
            port=port,
            database=database,
        )
        engine = create_engine(conn_url)
        return engine

    # ── PostgreSQL ────────────────────────────────────────────────────────────
    elif dbtype_lower in ('postgresql', 'postgres', 'pg'):
        if not (uid and pwd):
            raise ValueError("PostgreSQL requires uid and pwd.")
        port = port or 5432
        conn_url = URL.create(
            "postgresql+psycopg2",
            username=uid,
            password=pwd,
            host=server,
            port=port,
            database=database,
        )
        engine = create_engine(conn_url)
        return engine

    # ── SQLite ────────────────────────────────────────────────────────────────
    elif dbtype_lower in ('sqlite', 'sqlite3'):
        # `database` is the file path; server/uid/pwd are ignored
        engine = create_engine(f"sqlite:///{database}")
        return engine

    # ── Unsupported ───────────────────────────────────────────────────────────
    else:
        raise ValueError(
            f"Unsupported dbtype '{dbtype}'. "
            "Supported: mssql, sqlserver, sql_server, mariadb, mysql, postgresql, postgres, pg, sqlite"
        )  
    

def save_df(df, filename, file_format):
    out_path = os.path.join(args.data_store_path, filename)
    if file_format == 'parquet':
        df.to_parquet(out_path + '.parquet', index=False)
    else:
        df.to_csv(out_path + '.csv', index=False)
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Test SQLAlchemy DB connection.")
    parser.add_argument('--dbtype', type=str, required=True, help='Database type: mssql, postgresql, mysql, oracle')
    parser.add_argument('--server', type=str, required=True, help='Database server hostname or address')
    parser.add_argument('--database', type=str, required=True, help='Database name (for Oracle: service_name or SID)')
    parser.add_argument('--uid', type=str, help='Database user')
    parser.add_argument('--pwd', type=str, help='Database password')
    parser.add_argument('--port', type=int, help='Database port (optional)')
    parser.add_argument('--driver', type=str, help='ODBC driver for SQL Server (optional)')
    parser.add_argument('--sql_query_schema', type=str, default = 'dbo.', help='schema used alongside dtaabse in sql server. example: "dbo." ')
    parser.add_argument('--data_store_path', type=str, default = './data_folder/', help='Path in which the data files get stored. Default path: %(deafult)s ')
    parser.add_argument("--file_name_prefix", type =str, default ='omop', help ='saved filenames prefix. Default: %(default)s')
    parser.add_argument( '--file_format',  type=str, default='csv',  choices=['csv', 'parquet'], 
    help='Output file format for data exports (csv or parquet). Default: csv'
    )

 
    args = parser.parse_args()

    # engine = get_alchemy_engine ("UMWWPRDOMOPDB01.ad.umassmed.edu", "{args.database}")

    ### Creating data folder
    folder_path = args.data_store_path
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    
    try:
        engine = get_alchemy_engine(
            dbtype=args.dbtype, 
            server=args.server, 
            database=args.database, 
            uid=args.uid, 
            pwd=args.pwd,
            port=args.port,
            driver=args.driver
        )
        with engine.connect() as conn:
            print("Connection established successfully.")
    except Exception as e:
        print(f"Connection failed: {e}")



    ############# SQL script cohort ########################
    cohort_inital_create_sql = f'''
    WITH LatestVisit AS (
        SELECT ---TOP 300
            person_id,
            visit_occurrence_id,
            visit_start_date AS index_date,
            ROW_NUMBER() OVER (PARTITION BY person_id ORDER BY visit_start_date DESC) AS rn
        FROM
        {args.database}.{args.sql_query_schema}visit_occurrence

        WHERE
            visit_start_date BETWEEN '2018-04-01' AND '2019-04-30'
    ),
    CohortInitial AS (
        SELECT
            lv.person_id,
            lv.visit_occurrence_id,
            lv.index_date,
            YEAR(lv.index_date) - p.year_of_birth AS age_at_index,
            0 AS outcome, -- Default to 0, will update later based on conditions,
            l.location_source_value,
            p.race_source_value as Race_CD,
            p.ethnicity_source_value as Hispanic_CD,
            c.concept_code as Sex_CD

        FROM
            LatestVisit lv
            JOIN {args.database}.{args.sql_query_schema}person p ON lv.person_id = p.person_id
            LEFT JOIN {args.database}.{args.sql_query_schema}location l
                ON p.location_id = l.location_id
            LEFT JOIN {args.database}.{args.sql_query_schema}concept c on c.concept_id = p.gender_concept_id
        WHERE
            lv.rn = 1
            AND YEAR(lv.index_date) - p.year_of_birth >= 18
    ),
    PositiveOutcome AS (
        SELECT
            person_id,
            MIN(event_date) AS outcome_date
        FROM (
            SELECT
                co.person_id,
                MIN(co.condition_start_date) AS event_date
            FROM
                {args.database}.{args.sql_query_schema}condition_occurrence co
            INNER JOIN  {args.database}.{args.sql_query_schema}concept cc ON co.condition_concept_id = cc.concept_id
            WHERE
                cc.concept_name = 'Type 2 diabetes mellitus'
            GROUP BY
                co.person_id
            HAVING
                COUNT(DISTINCT co.condition_start_date) >= 2
            UNION ALL
            SELECT
                me.person_id,
                MIN(me.measurement_date) AS event_date
            FROM
                {args.database}.{args.sql_query_schema}measurement me
            INNER JOIN  {args.database}.{args.sql_query_schema}concept mc ON me.measurement_concept_id = mc.concept_id
            WHERE
                mc.concept_name = 'Hemoglobin A1c/Hemoglobin.total in Blood'
                AND me.value_as_number >= 6.5
            GROUP BY
                me.person_id
            HAVING
                COUNT(DISTINCT me.measurement_date) >= 2
        ) sub
        GROUP BY
            person_id
    ),
    ExclusionCriteria AS (
        SELECT
            ci.person_id
        FROM
            CohortInitial ci
            LEFT JOIN PositiveOutcome po ON ci.person_id = po.person_id
            LEFT JOIN death d ON ci.person_id = d.person_id
        WHERE
            (po.outcome_date IS NOT NULL AND po.outcome_date <= ci.index_date)
            OR (po.outcome_date IS NOT NULL AND d.death_date IS NOT NULL AND d.death_date < po.outcome_date)
    ),
    FinalCohort AS (
        SELECT
            ci.person_id,
            ci.visit_occurrence_id,
            ci.index_date,
            ci.age_at_index as Age_at_index,
            
            ci.location_source_value as Zipcode,
            ci.Race_CD,
            ci.Hispanic_CD,
            ci.Sex_CD,
            d.death_date,
            CASE
                WHEN po.outcome_date IS NOT NULL AND po.outcome_date <= DATEADD(YEAR, 5, ci.index_date)
                    THEN 1 ELSE 0
            END AS Outcome
        FROM
            CohortInitial ci
            LEFT JOIN PositiveOutcome po ON ci.person_id = po.person_id
            LEFT JOIN death d ON ci.person_id = d.person_id
        WHERE
            ci.person_id NOT IN (SELECT person_id FROM ExclusionCriteria)
    )
        SELECT * INTO {args.database}.{args.sql_query_schema}FinalCohort; --Writing table into the database
    '''
    ### above we are writing inital cohort to database (write function)

    print("Creating initial cohort file")
    try:
        with engine.connect() as conn:
            conn.execute(text(cohort_inital_create_sql))
            print("Query executed successfully.")
    except Exception as e:
        print(f'Error: {e}')
   
    print("Starting data retrieval")
    cohort_sql = f'''
    select * from {args.database}.{args.sql_query_schema}FinalCohort;
    '''


    df = pd.read_sql(cohort_sql, con=engine)
    # df.to_csv('{args.data_store_path}_patient_details.csv', index=False)
    save_df(df, args.file_name_prefix+ '_patient_details', args.file_format)
    print('--- Cohort creation completed---')

    patient_ids = df['person_id'].apply(str).tolist()
    pids_placeholder = ', '.join(patient_ids)

    ### medications ###
    med_sql = f'''
    SELECT
    x.person_id, 
    x.visit_occurrence_id, 
    x.drug_source_concept_id,
    a.ancestor_concept_id, 
    acon.concept_class_id AS ancestor_concept_class_id,
    a.min_levels_of_separation, 
    a.max_levels_of_separation,
    acon.concept_name AS ancestor_drug_concept_name,
    con.concept_code,
    con.concept_name,
    x.drug_exposure_start_date
FROM 
    {args.database}.{args.sql_query_schema}drug_exposure x
    JOIN {args.database}.{args.sql_query_schema}FinalCohort fc
        ON x.person_id = fc.person_id
    JOIN {args.database}.{args.sql_query_schema}CONCEPT_ANCESTOR a
        ON x.drug_source_concept_id = a.descendant_concept_id
    JOIN {args.database}.{args.sql_query_schema}CONCEPT acon
        ON a.ancestor_concept_id = acon.concept_id
    JOIN {args.database}.{args.sql_query_schema}CONCEPT con
        ON x.drug_source_concept_id = con.concept_id
WHERE 
    x.drug_exposure_start_date >= DATEADD(YEAR, -1, fc.index_date)
    AND x.drug_exposure_start_date <= fc.index_date
    AND acon.vocabulary_id = 'RxNorm'
    AND acon.standard_concept = 'S'
    AND acon.concept_class_id = 'Ingredient'
    AND a.max_levels_of_separation > 0
    ;
    '''
    meds_df = pd.read_sql(med_sql,con = engine)
    # meds_df.to_csv("{args.data_store_path}omop_medications.csv",index = False)
    save_df(meds_df, args.file_name_prefix+ '_medications', args.file_format)
    print('--- Retrieving medications data completed---')
    


    ### lab results ####
    lab_sql = f'''
    select  m.person_id,
        m.measurement_date,
        m.measurement_source_value, -- loinc code
        m.value_as_number,
        m.unit_source_value,
        c.*
    from 
        {args.database}.{args.sql_query_schema}measurement m
    join  {args.database}.{args.sql_query_schema}concept c 
        on m.measurement_concept_id = c.concept_id
    JOIN {args.database}.{args.sql_query_schema}FinalCohort fc
        ON m.person_id = fc.person_id

    WHERE c.vocabulary_id = 'LOINC'
    AND (m.measurement_date >= DATEADD(YEAR, -1, fc.index_date)
    AND m.measurement_date <= fc.index_date)
    '''
    lab_df = pd.read_sql(lab_sql,con = engine)
    save_df(lab_df, args.file_name_prefix+ '_labresults', args.file_format)
    print('--- Retrieving lab results data completed---')


    observation_vitals_sql = f'''
    select  ov.person_id,
    ov.observation_date,
    ov.observation_source_value, -- loinc code
    ov.value_as_number
    
    from 
        {args.database}.{args.sql_query_schema}Observation_Vitals ov
    JOIN {args.database}.{args.sql_query_schema}FinalCohort fc
        ON ov.person_id = fc.person_id
    where person_id IN ({pids_placeholder}) 
    and value_as_string NOT LIKE '%/%'
    and (ov.observation_date >= DATEADD(YEAR, -1, fc.index_date)
    AND ov.observation_date <= fc.index_date)
    '''
    lab_vitals_df = pd.read_sql(observation_vitals_sql,con = engine)
    save_df(lab_vitals_df, args.file_name_prefix+ '_labresults_vitals', args.file_format)
    print('--- Retrieving lab results vitals data completed---')



    #### diagnoses codes ####
    dx_sql = f'''
    select  m.person_id,
        m.condition_start_date,
        c.concept_code, -- ICD10 or ICD9 code
        c.vocabulary_id
    from {args.database}.{args.sql_query_schema}condition_occurrence m
    join {args.database}.{args.sql_query_schema}concept c on m.condition_source_concept_id = c.concept_id
    JOIN {args.database}.{args.sql_query_schema}FinalCohort fc
        ON m.person_id = fc.person_id

    WHERE c.vocabulary_id = 'LOINC'
    AND (m.measurement_date >= DATEADD(YEAR, -1, fc.index_date)
    AND m.measurement_date <= fc.index_date)
    '''
    dx_df = pd.read_sql(dx_sql,con = engine)
    save_df(lab_vitals_df, args.file_name_prefix+ '_diagnoses', args.file_format)

    print('--- Retrieving diagnoses data completed---')