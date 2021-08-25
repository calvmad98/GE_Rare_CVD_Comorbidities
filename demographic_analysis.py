import labkey
import pandas as pd
import numpy as np

# CARRIER AND NON-CARRIER DEMOGRAPHIC ANALYSIS

# Loading in table rare_diseases_participant_disease (info on disease onset and type)

labkey_server = "labkey-embassy.gel.zone"  
project_name = "main-programme/main-programme_v12_2021-05-06"  
context_path = "labkey"   
schema_name = "lists"  
query_name = "rare_diseases_participant_disease"  

# Create an object that will let us connect to the LabKey databases
server_context = labkey.utils.create_server_context(
    labkey_server, project_name, context_path, use_ssl=True
)
disease_results = labkey.query.select_rows(server_context, schema_name, query_name, max_rows=200000)
disease_results = pd.DataFrame(disease_results["rows"])
disease_results.head()

# Loading table with height and weight measurements
labkey_server = "labkey-embassy.gel.zone"  
project_name = "main-programme/main-programme_v12_2021-05-06"  
context_path = "labkey"   
schema_name = "lists"  
query_name = "rare_diseases_gen_measurement"  

# Create an object that will let us connect to the LabKey databases
server_context = labkey.utils.create_server_context(
    labkey_server, project_name, context_path, use_ssl=True
)
results = labkey.query.select_rows(server_context, schema_name, query_name, max_rows=200000)
gen_measure = pd.DataFrame(results["rows"])

# Loading in table cancer_risk_factor_general (info on participant height and weight)
# Specify what we are connecting to, and what schema and tables we want
labkey_server = "labkey-embassy.gel.zone"  
project_name = "main-programme/main-programme_v12_2021-05-06"  
context_path = "labkey"   
schema_name = "lists"  
query_name = "cancer_risk_factor_general"  

# Create an object that will let us connect to the LabKey databases
server_context = labkey.utils.create_server_context(
    labkey_server, project_name, context_path, use_ssl=True
)
cancer_measure = labkey.query.select_rows(server_context, schema_name, query_name, max_rows=200000)
cancer_measure = pd.DataFrame(cancer_measure["rows"])
cancer_measure.head()

# Reading in carriers data
carriers = pd.read_csv("all-gene_carrier.txt", sep = "\t")
carriers.head()

# Subsetting genes of interest
my_carriers = carriers.loc[(carriers['gene'] == 'MYBPC3') | 
                                     (carriers['gene'] == 'MYH7') |
                                     (carriers['gene'] == 'TNNT2') |
                                     (carriers['gene'] == 'TNNI3') |
                                     (carriers['gene'] == 'TPM1') |
                                     (carriers['gene'] == 'MYL3') |
                                     (carriers['gene'] == 'ACTC1') |
                                     (carriers['gene'] == 'PRKAG2') |
                                     (carriers['gene'] == 'GLA') |
                                     (carriers['gene'] == 'MYL2') |
                                     (carriers['gene'] == 'LMNA') |
                           (carriers['gene'] == 'PKP2') |
                           (carriers['gene'] == 'DSP') |
                           (carriers['gene'] == 'DSC2') |
                           (carriers['gene'] == 'TMEM43') |
                           (carriers['gene'] == 'DSG2') |
                           (carriers['gene'] == 'KCNQ1') |
                           (carriers['gene'] == 'KCNH2') |
                           (carriers['gene'] == 'SCN5A')]

my_carriers = my_carriers.sort_values('participant_id')
my_carriers = my_carriers.reset_index()

# Reading in all participants with genomic information
participants = pd.read_csv("participants_in_gvcf.txt", sep = "\t")

# Finding participants in participants file with the genes of interest
my_carriers = participants.loc[participants['participant_id'].isin(my_carriers.participant_id)] 
my_carriers.to_csv("all_carrier_demographics.txt", index = False)
print("Data written to file!")

# Looking at carriers with genes of interest
gene_carriers = pd.read_csv("all-gene_carrier.txt", sep = "\t")

# Creating died column indicating if participant has a date of death
my_carriers['died'] = np.where(pd.notnull(my_carriers['date_of_death']), 1, None)

# Extracting the year of death from the date of death column
my_carriers['year_of_death'] = pd.DatetimeIndex(my_carriers['date_of_death']).year

# Calculate age from either birth date or death date
my_carriers['age'] = ((np.where(pd.isnull(my_carriers['date_of_death']), 2021, my_carriers['year_of_death']) - my_carriers['year_of_birth']))

my_carriers = my_carriers

# Observing if any duplicates in demographic table
duplicates = my_carriers[my_carriers.duplicated(['participant_id'], keep = False)]
len(duplicates.participant_id.unique())

# Extracting height and weight data, cleaning
gen_rd = gen_measure.loc[gen_measure['participant_id'].isin(my_carriers.participant_id.unique())] 
height_rd = gen_rd.loc[(gen_rd['measurement_classification'] == 'Height')]
weight_rd = gen_rd.loc[(gen_rd['measurement_classification'] == 'Weight')]
height_rd.loc[(height_rd.measurement > 500),'measurement']= None
height_rd.loc[height_rd['measurement'] > 5, 'measurement'] = height_rd['measurement'].div(100)
height_rd.loc[(height_rd.measurement == 0.00),'measurement']= None
weight_rd.loc[(weight_rd.measurement == 0.00),'measurement']= None

# Finding height and weight data for carriers that are cancer participants
gen_c = cancer_measure.loc[cancer_measure['participant_id'].isin(my_carriers.participant_id.unique())] 
cancer_h_w = gen_c[['participant_id', 'height', 'weight']]
cancer_h_w.loc[(cancer_h_w.height > 500),'height']= None
cancer_h_w.loc[cancer_h_w['height'] > 5, 'height'] = cancer_h_w['height'].div(100)

# Cleaning height and weight columns
height_rd = height_rd.rename(columns = {'measurement':'height_rare_disease'})
weight_rd = weight_rd.rename(columns = {'measurement':'weight_rare_disease'})
cancer_h_w = cancer_h_w.rename(columns = {'height':'height_cancer',
                                         'weight':'weight_cancer'})

# Information on age of onset from disease_results table 
disease_dems = disease_results.loc[disease_results['participant_id'].isin(my_carriers.participant_id.unique())] 

# Subsetting columns of interest
carrier_demographics = my_carriers[['participant_id', 
                        'participant_type', 
                        'participant_phenotypic_sex', 
                        'age', 
                        'participant_ethnic_category_simple', 
                        'ancestry',
                        'died',
                        'date_of_death', 
                        'imd_newest', 
                        'imd_lowest', 
                        'genome_build', 
                        'related', 
                        'pred_african_ancestries',
                        'pred_south_asian_ancestries', 
                        'pred_east_asian_ancestries',
                        'pred_european_ancestries', 
                        'pred_american_ancestries']] 

# Merging all wanted information
# Information on age of onset and recruitment disease
carrier_demographics = pd.merge(carrier_demographics, disease_dems[['participant_id', 
                                                                    'age_of_onset', 
                                                                    'specific_disease']], 
                            on='participant_id', how='left' )

# Height info
carrier_demographics = pd.merge(carrier_demographics, height_rd[['participant_id', 'height_rare_disease']], 
                            on='participant_id', how='left' )
# Weight info
carrier_demographics = pd.merge(carrier_demographics, weight_rd[['participant_id', 'weight_rare_disease']], 
                            on='participant_id', how='left' )
# Cancer height and weight info
carrier_demographics = pd.merge(carrier_demographics, cancer_h_w[['participant_id', 'height_cancer', 
                                                                        'weight_cancer']], 
                                   on='participant_id', how='left' )

# Cleaning height and weight columns
carrier_demographics['height'] = np.where(pd.isnull(carrier_demographics['height_cancer'].values), 
                                          carrier_demographics['height_rare_disease'], 
                                          carrier_demographics['height_cancer'])
carrier_demographics['weight'] = np.where(pd.isnull(carrier_demographics['weight_cancer'].values), 
                                          carrier_demographics['weight_rare_disease'], 
                                          carrier_demographics['weight_cancer'])

# Dropping original columns for height and weight
carrier_demographics = carrier_demographics.drop(columns = ['height_rare_disease', 'height_cancer', 
                                                            'weight_rare_disease', 'weight_cancer'])

# Merging information on pathogenic variant
carrier_demographics = pd.merge(carrier_demographics, gene_carriers[['participant_id', 
                                                       'gene', 
                                                       'ID', 
                                                       'genotype', 
                                                       'CLNSIG', 
                                                       'CLNVC', 
                                                       'MC', 
                                                       'source']], 
                            on='participant_id', how='left' )

# Loading in table cancer_analysis (info on specific cancer)
labkey_server = "labkey-embassy.gel.zone" 
project_name = "main-programme/main-programme_v12_2021-05-06"  
context_path = "labkey"  
schema_name = "lists"  
query_name = "cancer_analysis"  
# Create an object that will let us connect to the LabKey databases. This does not change.
server_context = labkey.utils.create_server_context(
    labkey_server, project_name, context_path, use_ssl=True
)
cancer_disease = labkey.query.select_rows(server_context, schema_name, query_name, max_rows=200000)
cancer_disease = pd.DataFrame(cancer_disease["rows"])

# Adding specific disease information for cancer carriers
m = cancer_disease.set_index('participant_id')['disease_type'].to_dict()
v = carrier_demographics.filter(items='specific_disease')
carrier_demographics[v.columns] = v.replace(m)

# Creating bins for age and age of onset information
carrier_demographics['age_category'] = pd.cut(x=carrier_demographics['age'], bins=[0,25,50,75,100,115])
carrier_demographics['onset_category'] = pd.cut(x=carrier_demographics['age_of_onset'], bins=[-1000,0,10,20,30,40,50,115])

# Sorting values by participant ID
carrier_demographics = carrier_demographics.sort_values('participant_id')

# Observing if any duplicates in demographic table
duplicates = carrier_demographics[carrier_demographics.duplicated(keep = False)]
len(duplicates.participant_id.unique())

# Dropping duplicate rows
carrier_demographics = carrier_demographics.drop_duplicates()

# Writing to file
carrier_demographics.to_csv("full_carriers_info.txt", index = False)
print("Data written to file!")

# Writing carrier IDs to file
carrier_ids = carrier_demographics[['participant_id']]
carrier_ids = carrier_ids.drop_duplicates()
carrier_ids.to_csv("carrier_ids.txt", index = False)
print("Data written to file!")

# # Splitting into Diseases
# Creating a carriers dataframe for each disease of interest 

HCM_carriers = carrier_demographics.loc[(carrier_demographics['gene'] == 'MYBPC3') | 
    (carrier_demographics['gene'] == 'MYH7') |
    (carrier_demographics['gene'] == 'TNNI3') |
    (carrier_demographics['gene'] == 'TPM1') |
    (carrier_demographics['gene'] == 'MYL3') |
    (carrier_demographics['gene'] == 'ACTC1') |
    (carrier_demographics['gene'] == 'MYL2') |
    (carrier_demographics['gene'] == 'PRKAG2')]
DCM_carriers = carrier_demographics.loc[(carrier_demographics['gene'] == 'TNNT2') | 
                            (carrier_demographics['gene'] == 'LMNA')]
FD_carriers = carrier_demographics.loc[carrier_demographics['gene'] == 'GLA']
ARVC_carriers = carrier_demographics.loc[(carrier_demographics['gene'] == 'PKP2') |
                             (carrier_demographics['gene'] == 'DSP') |
                             (carrier_demographics['gene'] == 'DSC2') |
                             (carrier_demographics['gene'] == 'TMEM43') |
                             (carrier_demographics['gene'] == 'DSG2')]
LQTS_BS_carriers = carrier_demographics.loc[(carrier_demographics['gene'] == 'KCNQ1') |
                                (carrier_demographics['gene'] == 'KCNH2') |
                                (carrier_demographics['gene'] == 'SCN5A')]

HCM_carriers = HCM_carriers.sort_values('participant_id')
DCM_carriers = DCM_carriers.sort_values('participant_id')
FD_carriers = FD_carriers.sort_values('participant_id')
ARVC_carriers = ARVC_carriers.sort_values('participant_id')
LQTS_BS_carriers = LQTS_BS_carriers.sort_values('participant_id')

# Dictionary containing all disease dataframes
carrier_list = {'Hypertrophic Cardiomyopathy': HCM_carriers, 
              'Dilated Cardiomyopathy': DCM_carriers, 
              'Fabry Disease': FD_carriers, 
              'Arrhythmogenic Right Ventricular Cardiomyopathy': ARVC_carriers, 
              'Long QT / Brugada Syndrome': LQTS_BS_carriers}

# Looking at duplicates in dataframes
for key in carrier_list:
    duplicates = carrier_list[key][carrier_list[key].duplicated()]
    print(key, 'duplicates: ', len(duplicates))
    print('----------------------------------------------')

# Ensuring all dataframes have one row per participant, and removing participants with missing sex
HCM_carriers = HCM_carriers.drop_duplicates(['participant_id'])
HCM_carriers = HCM_carriers.loc[HCM_carriers['participant_phenotypic_sex'] != 'Indeterminate']
DCM_carriers = DCM_carriers.drop_duplicates(['participant_id'])
DCM_carriers = DCM_carriers.loc[DCM_carriers['participant_phenotypic_sex'] != 'Indeterminate']
FD_carriers = FD_carriers.drop_duplicates(['participant_id'])
FD_carriers = FD_carriers.loc[FD_carriers['participant_phenotypic_sex'] != 'Indeterminate']
ARVC_carriers = ARVC_carriers.drop_duplicates(['participant_id'])
ARVC_carriers = ARVC_carriers.loc[ARVC_carriers['participant_phenotypic_sex'] != 'Indeterminate']
LQTS_BS_carriers = LQTS_BS_carriers.drop_duplicates(['participant_id'])
LQTS_BS_carriers = LQTS_BS_carriers.loc[LQTS_BS_carriers['participant_phenotypic_sex'] != 'Indeterminate']

# Dictionary containing all disease dataframes
carrier_list = {'Hypertrophic Cardiomyopathy': HCM_carriers, 
              'Dilated Cardiomyopathy': DCM_carriers, 
              'Fabry Disease': FD_carriers, 
              'Arrhythmogenic Right Ventricular Cardiomyopathy': ARVC_carriers, 
              'Long QT / Brugada Syndrome': LQTS_BS_carriers}

# Writing to files
HCM_carriers.to_csv("HCM_carrier_demographics.txt", index = False)
DCM_carriers.to_csv("DCM_carrier_demographics.txt", index = False)
FD_carriers.to_csv("FD_carrier_demographics.txt", index = False)
ARVC_carriers.to_csv("ARVC_carrier_demographics.txt", index = False)
LQTS_BS_carriers.to_csv("LQTS_BS_carrier_demographics.txt", index = False)
print("Data written to file!")

# Creating function to count each gender 
def findGenderCount(carrier_dataframe, id_list):
    
    male_df = carrier_dataframe.loc[(carrier_dataframe['participant_id'].isin(id_list)) & 
                      (carrier_dataframe['participant_phenotypic_sex'] == 'Male')]
    female_df = carrier_dataframe.loc[(carrier_dataframe['participant_id'].isin(id_list)) & 
                      (carrier_dataframe['participant_phenotypic_sex'] == 'Female')]
    males = len(male_df.participant_id.unique()) 
    females = len(female_df.participant_id.unique())
    
    final_table = pd.DataFrame({'Males' : males,
                               'Females' : females},
                              index = ['Total Count']) 
    return (final_table);
for key in carrier_list:
    print(key)
    print(findGenderCount(carrier_demographics, carrier_list[key].participant_id))
    print('------------------------------------')

# Creating function to count each ethnicity 
def findEthnicity(carrier_dataframe, id_list):
    temp_df = carrier_dataframe.loc[(carrier_dataframe['participant_id'].isin(id_list))]
    temp_df = temp_df.drop_duplicates(['participant_id'])
    value_dataframe = temp_df['participant_ethnic_category_simple'].value_counts().rename_axis('Ethnicity').reset_index(name='count')
    return value_dataframe;
for key in carrier_list:
    print(key)
    print(findEthnicity(carrier_demographics, carrier_list[key].participant_id))
    print('-----------------------------------')

# Creating function to count each age category 
def findAgeCategory(carrier_dataframe, id_list):
    temp_df = carrier_dataframe.loc[(carrier_dataframe['participant_id'].isin(id_list))]
    temp_df = temp_df.drop_duplicates(['participant_id'])
    value_df = temp_df['age_category'].value_counts().rename_axis('Age Category').reset_index(name='count')
    return value_df;
for key in carrier_list:
    print(key)
    print(findAgeCategory(carrier_demographics, carrier_list[key].participant_id))
    print('------------------------------')

# Creating function to count each ancestry 
def findAncestry(carrier_dataframe, id_list):
    temp_df = carrier_dataframe.loc[(carrier_dataframe['participant_id'].isin(id_list))]
    temp_df = temp_df.drop_duplicates(['participant_id'])
    value_df = temp_df['ancestry'].value_counts().rename_axis('Ancestry').reset_index(name='count')
    return value_df;
                                                                                                                                                                                                               
for key in carrier_list:
    print(key)
    print(findAncestry(carrier_demographics, carrier_list[key].participant_id))
    print('-------------------------------')

# Creating function to find median height
def findMedianHeight(carrier_dataframe, id_list):
    temp_df = carrier_dataframe.loc[(carrier_dataframe['participant_id'].isin(id_list))]
    temp_df = temp_df.drop_duplicates(['participant_id'])
    return temp_df['height'].median();
for key in carrier_list:
    print(key)
    print(findMedianHeight(carrier_demographics, carrier_list[key].participant_id))
    print('----------------------')

# Creating function to find median weight
def findMedianWeight(carrier_dataframe, id_list):
    temp_df = carrier_dataframe.loc[(carrier_dataframe['participant_id'].isin(id_list))]
    temp_df = temp_df.drop_duplicates(['participant_id'])
    return temp_df['weight'].median();
for key in carrier_list:
    print(key)
    print(findMedianWeight(carrier_demographics, carrier_list[key].participant_id))
    print('----------------------')

# Creating function to count age of onset bins
def findAgeOfOnset(carrier_dataframe, id_list):
    temp_df = carrier_dataframe.loc[(carrier_dataframe['participant_id'].isin(id_list))]
    temp_df = temp_df.drop_duplicates(['participant_id'])
    value_df = temp_df['onset_category'].value_counts().rename_axis('Age of Onset').reset_index(name='count')
    return value_df;
for key in carrier_list:
    print(key)
    print(findAgeOfOnset(carrier_demographics, carrier_list[key].participant_id))
    print('---------------------------------')


# # Building Demographic Table
# Creating dataframe with demographic information
# Total carriers per disease

carrier_demographics_table = pd.DataFrame([[len(HCM_carriers.participant_id.unique()),
                                            len(DCM_carriers.participant_id.unique()), 
                                            len(FD_carriers.participant_id.unique()),
                                            len(ARVC_carriers.participant_id.unique()), 
                                            len(LQTS_BS_carriers.participant_id.unique())]],
                                         index= ['Total Probands and Relatives'],
                                         columns = ['HCM', 'DCM', 'FD', 'ARVC', 'LQTS_BS'])

# Adding gender information
gender_table = pd.DataFrame({'HCM': [125, 136],
                             'DCM': [78, 77],
                             'FD': [9, 16],
                             'ARVC': [97, 155],
                             'LQTS_BS': [162, 195]}, 
                            index=pd.MultiIndex.from_tuples(
                                [('Gender','Males'), ('Gender','Females')], 
                                names=['Demographic', 'Subgroup']))
carrier_demographics_table = pd.concat([carrier_demographics_table, gender_table])

# Adding ethnicity information
ethnicity_table = pd.DataFrame({'HCM': [175, 46, 17, 11, 9, '<5'],
                                'DCM': [73, 28, 9, 26, 6, 13],
                                'FD': [12, 5, 8, 0, 0, 0],
                             'ARVC': [180, 36, 15, 11, '<5', 8],
                             'LQTS_BS': [231, 55, 38, 17, 10, 6]}, 
                            index=pd.MultiIndex.from_tuples(
                                [('Ethnicity','White'), ('Ethnicity','Not known / Missing'), 
                                 ('Ethnicity','Asian or Asian British'), ('Ethnicity','Black or Black British'),
                                 ('Ethnicity','Other'), ('Ethnicity','Mixed')], 
                                names=['Demographic', 'Subgroup']))
carrier_demographics_table = pd.concat([carrier_demographics_table, ethnicity_table])

# Adding ancestry information
ancestry_table = pd.DataFrame({'HCM': [206, 55],
                               'DCM': [91, 64],
                               'FD': [16, 9],
                             'ARVC': [199, 53],
                             'LQTS_BS': [251, 106]}, 
                            index=pd.MultiIndex.from_tuples(
                                [('Ancestry','European'), ('Ancestry','Other')], 
                                names=['Demographic', 'Subgroup']))
carrier_demographics_table = pd.concat([carrier_demographics_table, ancestry_table])

# Adding death information
died_table = pd.DataFrame({'HCM': [12],
                           'DCM': [8],
                           'FD': ['<5'],
                             'ARVC': [14],
                             'LQTS_BS': [23]}, 
                            index=pd.MultiIndex.from_tuples(
                                [('Death','Count')], 
                                names=['Demographic', 'Subgroup']))
carrier_demographics_table = pd.concat([carrier_demographics_table, died_table])

# Adding height information
height_table = pd.DataFrame({'HCM': [1.62],
                             'DCM': [1.15],
                             'FD': [0.59],
                             'ARVC': [1.68],
                             'LQTS_BS': [1.68]}, 
                            index=pd.MultiIndex.from_tuples(
                                [('Height','Median (m)')], 
                                names=['Demographic', 'Subgroup']))
carrier_demographics_table = pd.concat([carrier_demographics_table, height_table])

# Adding weight information
weight_table = pd.DataFrame({'HCM': [70.3],
                             'DCM': [21.9],
                             'FD': [6.52],
                             'ARVC': [61.9],
                             'LQTS_BS': [75.8]}, 
                            index=pd.MultiIndex.from_tuples(
                                [('Weight','Median (kg)')], 
                                names=['Demographic', 'Subgroup']))
carrier_demographics_table = pd.concat([carrier_demographics_table, weight_table])

# Adding age information
age_table = pd.DataFrame({'HCM': [42, 107, 97, 15, 0],
                          'DCM': [42, 67, 36, 9, 0],
                          'FD': ['<5', 7, 11, '<5', 0],
                             'ARVC': [42, 100, 94, 15, 0],
                             'LQTS_BS': [73, 139, 116, 28, 0]}, 
                            index=pd.MultiIndex.from_tuples(
                                [('Age','0-25'), ('Age','25-50'), 
                                 ('Age','50-75'), ('Age','75-100'),
                                 ('Age','100-115')], 
                                names=['Demographic', 'Subgroup']))
carrier_demographics_table = pd.concat([carrier_demographics_table, age_table])

# Adding age of onset information
onset_table = pd.DataFrame({'HCM': [45, 16, 17, 18, 24, 17, 16],
                            'DCM': [38, 11, 8, 8, 5, '<5', 7],
                            'FD': [6, 0, 0, '<5', '<5', 0, '<5'],
                             'ARVC': [48, 14, 9, 13, 10, 11, 5],
                             'LQTS_BS': [56, 25, 16, 14, 10, 15, 10]}, 
                            index=pd.MultiIndex.from_tuples(
                                [('Age of Onset','Prenatal'), ('Age of Onset','0-10'), ('Age of Onset','10-20'), 
                                 ('Age of Onset','20-30'), ('Age of Onset','30-40'),
                                 ('Age of Onset','40-50'), ('Age of Onset', '50-115')], 
                                names=['Demographic', 'Subgroup']))
carrier_demographics_table = pd.concat([carrier_demographics_table, onset_table])

# Renaming columns
carrier_demographics_table = carrier_demographics_table.rename(columns = {'HCM':'HCM_carriers',
                                                                          'DCM':'DCM_carriers',
                                                                          'FD':'FD_carriers',
                                                                          'ARVC':'ARVC_carriers',
                                                                          'LQTS_BS':'LQTS_BS_carriers'})

# Creating category column to indicate what each row represents
carrier_demographics_table['category'] = ['Total Probands and Relatives' , '(Gender, Males)','(Gender, Females)',
                                            '(Ethnicity, White)','(Ethnicity, Not known / Missing)',
                                            '(Ethnicity, Asian or Asian British)','(Ethnicity, Black or Black British)',
                                            '(Ethnicity, Other)','(Ethnicity, Mixed)','(Ancestry, European)',
                                            '(Ancestry, Other)','(Death, Count)','(Height, Median (m))',
                                            '(Weight, Median (kg))','(Age, 0-25)','(Age, 25-50)','(Age, 50-75)',
                                            '(Age, 75-100)','(Age, 100-115)','(Age of Onset, Prenatal)',
                                            '(Age of Onset, 0-10)','(Age of Onset, 10-20)' ,'(Age of Onset, 20-30)',
                                            '(Age of Onset, 30-40)','(Age of Onset, 40-50)','(Age of Onset, 50-115)' ]

# Writing to file
carrier_demographics_table.to_csv("final_carrier_demographics_table.txt", index = False)
print("Data written to file!")

# NON-CARRIER DEMOGRAPHICS

# Loading in table rare_diseases_gen_measurement (info on participant height and weight)
labkey_server = "labkey-embassy.gel.zone"  
project_name = "main-programme/main-programme_v12_2021-05-06"
context_path = "labkey" 
schema_name = "lists" 
query_name = "rare_diseases_gen_measurement" 
server_context = labkey.utils.create_server_context(
    labkey_server, project_name, context_path, use_ssl=True
)
results = labkey.query.select_rows(server_context, schema_name, query_name, max_rows=200000)
gen_measure = pd.DataFrame(results["rows"])
gen_measure.head()

# Loading in table cancer_risk_factor_general (info on participant height and weight)
labkey_server = "labkey-embassy.gel.zone"  
project_name = "main-programme/main-programme_v12_2021-05-06"  
context_path = "labkey"  
schema_name = "lists"  
query_name = "cancer_risk_factor_general"  
server_context = labkey.utils.create_server_context(
    labkey_server, project_name, context_path, use_ssl=True
)
cancer_measure = labkey.query.select_rows(server_context, schema_name, query_name, max_rows=200000)
cancer_measure = pd.DataFrame(cancer_measure["rows"])
cancer_measure.head()

# Loading in table rare_diseases_participant_disease (info on age of onset)
labkey_server = "labkey-embassy.gel.zone"  
project_name = "main-programme/main-programme_v12_2021-05-06" 
context_path = "labkey"   
schema_name = "lists"  
query_name = "rare_diseases_participant_disease" 
server_context = labkey.utils.create_server_context(
    labkey_server, project_name, context_path, use_ssl=True
)
disease_results = labkey.query.select_rows(server_context, schema_name, query_name, max_rows=200000)
disease_results = pd.DataFrame(disease_results["rows"])
disease_results.head()

# Loading in information on each disease carriers
HCM_carriers = pd.read_csv("HCM_carriers.txt", sep = ",")
DCM_carriers = pd.read_csv("DCM_carriers.txt", sep = ",")
FD_carriers = pd.read_csv("FD_carriers.txt", sep = ",")
ARVC_carriers =  pd.read_csv("ARVC_carriers.txt", sep = ",")
LQTS_BS_carriers = pd.read_csv("LQTS_BS_carriers.txt", sep = ",")

# Reading in participants with genomics information
participants = pd.read_csv("participants_in_gvcf.txt", sep = "\t")
participants.head(5)

# Dictionary with carrier dataframes
carrier_list = {'Hypertrophic Cardiomyopathy': HCM_carriers, 
              'Dilated Cardiomyopathy': DCM_carriers, 
              'Fabry Disease': FD_carriers, 
              'Arrhythmogenic Right Ventricular Cardiomyopathy': ARVC_carriers, 
              'Long QT / Brugada Syndrome': LQTS_BS_carriers}

# Calculating total non-carriers per disease
print('Total number of participants: ', len(participants.participant_id.unique()))
print('')
for key in carrier_list:
    print('Total number of', key, 'non-carriers:', len(participants.participant_id.unique()) - len(carrier_list[key].participant_id.unique()))
    print('--------------------------------------------------')

# Identifying non-carriers for each disease
HCM_noncarriers = participants.loc[~participants['participant_id'].isin(HCM_carriers.participant_id.unique())]
DCM_noncarriers = participants.loc[~participants['participant_id'].isin(DCM_carriers.participant_id.unique())]
FD_noncarriers = participants.loc[~participants['participant_id'].isin(FD_carriers.participant_id.unique())]
ARVC_noncarriers = participants.loc[~participants['participant_id'].isin(ARVC_carriers.participant_id.unique())]
LQTS_BS_noncarriers = participants.loc[~participants['participant_id'].isin(LQTS_BS_carriers.participant_id.unique())]

# Non-carrier dataframes in dictionary
disease_lists = [HCM_noncarriers, DCM_noncarriers, FD_noncarriers, ARVC_noncarriers, LQTS_BS_noncarriers]

# Create column indicating if the participant has died or not
for i in range(len(disease_lists)):
     disease_lists[i]['died'] = np.where(pd.notnull(disease_lists[i]['date_of_death']), 1, None)

# Create Year of Death column
for i in range(len(disease_lists)):
    disease_lists[i]['year_of_death'] = pd.DatetimeIndex(disease_lists[i]['date_of_death']).year

# Calculate age
for i in range(len(disease_lists)):
    disease_lists[i]['age'] = ((np.where(pd.isnull(disease_lists[i]['date_of_death']), 2021, disease_lists[i]['year_of_death']) - disease_lists[i]['year_of_birth']))

for i in range(len(disease_lists)):
    duplicates = disease_lists[i][disease_lists[i].duplicated(['participant_id'], keep = False)]
    print(len(duplicates.participant_id.unique()))
    print('')

# Obtaining height and weight information
height_rd = gen_measure.loc[(gen_measure['measurement_classification'] == 'Height')]
weight_rd = gen_measure.loc[(gen_measure['measurement_classification'] == 'Weight')]
height_rd.loc[(height_rd.measurement == 0.00),'measurement']= None
weight_rd.loc[(weight_rd.measurement == 0.00),'measurement']= None

# Cleaning height info
height_rd.loc[(height_rd.measurement > 500),'measurement']= None
height_rd.loc[height_rd['measurement'] > 5, 'measurement'] = height_rd['measurement'].div(100)

# Height and weight info for cancer patients
cancer_h_w = cancer_measure[['participant_id', 'height', 'weight']]
cancer_h_w.loc[(cancer_h_w.height == 0.00),'height']= None
cancer_h_w.loc[(cancer_h_w.weight == 0.00),'weight']= None

# Cleaning cancer height info
cancer_h_w.loc[(cancer_h_w.height > 500),'height']= None
cancer_h_w.loc[cancer_h_w['height'] > 5, 'height'] = cancer_h_w['height'].div(100)

# Cleaning dataframes
height_rd = height_rd.rename(columns = {'measurement':'height_rare_disease'})
weight_rd = weight_rd.rename(columns = {'measurement':'weight_rare_disease'})
cancer_h_w = cancer_h_w.rename(columns = {'height':'height_cancer',
                                         'weight':'weight_cancer'})

# Merging height and weight info to disease dataframe
HCM_noncarriers = pd.merge(HCM_noncarriers, height_rd[['participant_id', 'height_rare_disease']], 
                            on='participant_id', how='left' )
HCM_noncarriers = pd.merge(HCM_noncarriers, weight_rd[['participant_id', 'weight_rare_disease']], 
                            on='participant_id', how='left' )
HCM_noncarriers = pd.merge(HCM_noncarriers, cancer_h_w[['participant_id', 'height_cancer', 'weight_cancer']], 
                            on='participant_id', how='left' )

DCM_noncarriers = pd.merge(DCM_noncarriers, height_rd[['participant_id', 'height_rare_disease']], 
                            on='participant_id', how='left' )
DCM_noncarriers = pd.merge(DCM_noncarriers, weight_rd[['participant_id', 'weight_rare_disease']], 
                            on='participant_id', how='left' )
DCM_noncarriers = pd.merge(DCM_noncarriers, cancer_h_w[['participant_id', 'height_cancer', 'weight_cancer']], 
                            on='participant_id', how='left' )

FD_noncarriers = pd.merge(FD_noncarriers, height_rd[['participant_id', 'height_rare_disease']], 
                            on='participant_id', how='left' )
FD_noncarriers = pd.merge(FD_noncarriers, weight_rd[['participant_id', 'weight_rare_disease']], 
                            on='participant_id', how='left' )
FD_noncarriers = pd.merge(FD_noncarriers, cancer_h_w[['participant_id', 'height_cancer', 'weight_cancer']], 
                            on='participant_id', how='left' )

ARVC_noncarriers = pd.merge(ARVC_noncarriers, height_rd[['participant_id', 'height_rare_disease']], 
                            on='participant_id', how='left' )
ARVC_noncarriers = pd.merge(ARVC_noncarriers, weight_rd[['participant_id', 'weight_rare_disease']], 
                            on='participant_id', how='left' )
ARVC_noncarriers = pd.merge(ARVC_noncarriers, cancer_h_w[['participant_id', 'height_cancer', 'weight_cancer']], 
                            on='participant_id', how='left' )

LQTS_BS_noncarriers = pd.merge(LQTS_BS_noncarriers, height_rd[['participant_id', 'height_rare_disease']], 
                            on='participant_id', how='left' )
LQTS_BS_noncarriers = pd.merge(LQTS_BS_noncarriers, weight_rd[['participant_id', 'weight_rare_disease']], 
                            on='participant_id', how='left' )
LQTS_BS_noncarriers = pd.merge(LQTS_BS_noncarriers, cancer_h_w[['participant_id', 'height_cancer', 'weight_cancer']], 
                            on='participant_id', how='left' )

# Recreating non-carrier dictionary
disease_lists = [HCM_noncarriers, DCM_noncarriers, FD_noncarriers, ARVC_noncarriers, LQTS_BS_noncarriers]

# Height Column as either rare disease or cancer height
for i in range(len(disease_lists)):
    disease_lists[i]['height'] = np.where(pd.isnull(disease_lists[i]['height_cancer'].values), 
                                          disease_lists[i]['height_rare_disease'], 
                                          disease_lists[i]['height_cancer'])

# Weight Column as either rare disease or cancer height
for i in range(len(disease_lists)):
    disease_lists[i]['weight'] = np.where(pd.isnull(disease_lists[i]['weight_cancer'].values), 
                                          disease_lists[i]['weight_rare_disease'], 
                                          disease_lists[i]['weight_cancer'])

# Dropping original height and weight columns
HCM_noncarriers = HCM_noncarriers.drop(columns = ['height_rare_disease', 'height_cancer', 
                                                        'weight_rare_disease', 'weight_cancer'])
DCM_noncarriers = DCM_noncarriers.drop(columns = ['height_rare_disease', 'height_cancer', 
                                                        'weight_rare_disease', 'weight_cancer'])
FD_noncarriers = FD_noncarriers.drop(columns = ['height_rare_disease', 'height_cancer', 
                                                        'weight_rare_disease', 'weight_cancer'])
ARVC_noncarriers = ARVC_noncarriers.drop(columns = ['height_rare_disease', 'height_cancer', 
                                                        'weight_rare_disease', 'weight_cancer'])
LQTS_BS_noncarriers = LQTS_BS_noncarriers.drop(columns = ['height_rare_disease', 'height_cancer', 
                                                        'weight_rare_disease', 'weight_cancer'])

# Recreating non-carrier dictionary
disease_lists = [HCM_noncarriers, DCM_noncarriers, FD_noncarriers, ARVC_noncarriers, LQTS_BS_noncarriers]

# Merging Age of Onset information
HCM_noncarriers = pd.merge(HCM_noncarriers, disease_results[['participant_id', 
                                                             'age_of_onset',
                                                             'specific_disease']], 
                                on='participant_id', how='left' )
DCM_noncarriers = pd.merge(DCM_noncarriers, disease_results[['participant_id', 
                                                             'age_of_onset',
                                                             'specific_disease']], 
                                on='participant_id', how='left' )
FD_noncarriers = pd.merge(FD_noncarriers, disease_results[['participant_id', 
                                                             'age_of_onset',
                                                             'specific_disease']], 
                                on='participant_id', how='left' )
ARVC_noncarriers = pd.merge(ARVC_noncarriers, disease_results[['participant_id', 
                                                             'age_of_onset',
                                                             'specific_disease']], 
                                on='participant_id', how='left' )
LQTS_BS_noncarriers = pd.merge(LQTS_BS_noncarriers, disease_results[['participant_id', 
                                                             'age_of_onset',
                                                             'specific_disease']], 
                                on='participant_id', how='left' )

# Loading in table cancer_analysis (info on specific disease for cancer participants)
labkey_server = "labkey-embassy.gel.zone" 
project_name = "main-programme/main-programme_v12_2021-05-06"
context_path = "labkey" 
schema_name = "lists" 
query_name = "cancer_analysis”
server_context = labkey.utils.create_server_context(
    labkey_server, project_name, context_path, use_ssl=True
)
cancer_disease = labkey.query.select_rows(server_context, schema_name, query_name, max_rows=200000)
cancer_disease = pd.DataFrame(cancer_disease["rows"])
cancer_disease.head()

# Recreating non-carrier dictionary
disease_lists = [HCM_noncarriers, DCM_noncarriers, FD_noncarriers, ARVC_noncarriers, LQTS_BS_noncarriers]

# Adding disease type to cancer patients
for i in range(len(disease_lists)): 
    m = cancer_disease.set_index('participant_id')['disease_type'].to_dict()
    v = disease_lists[i].filter(items='specific_disease')
    disease_lists[i][v.columns] = v.replace(m)

# Creating age bins
for i in range(len(disease_lists)): 
    disease_lists[i]['age_category'] = pd.cut(x=disease_lists[i]['age'], bins=[0,25,50,75,100,115])
    disease_lists[i]['onset_category'] = pd.cut(x=disease_lists[i]['age_of_onset'], bins=[-1000,0,10,20,30,40,50,115])

# Sorting dataframes by participant ID
for i in range(len(disease_lists)): 
    disease_lists[i] = disease_lists[i].sort_values('participant_id')

# Observing if any duplicates in demographic tables
for i in range(len(disease_lists)):
    duplicates = disease_lists[i][disease_lists[i].duplicated(keep = False)]
    print(len(disease_lists[i]))
    print(len(duplicates))
    print('')

# Dropping duplicate rows
for i in range(len(disease_lists)):
    disease_lists[i] = disease_lists[i].drop_duplicates()
    print(len(disease_lists[i]))

HCM_noncarriers = HCM_noncarriers.drop_duplicates(['participant_id'])
HCM_noncarriers = HCM_noncarriers.loc[HCM_noncarriers['participant_phenotypic_sex'] != 'Indeterminate']
DCM_noncarriers = DCM_noncarriers.drop_duplicates(['participant_id'])
DCM_noncarriers = DCM_noncarriers.loc[DCM_noncarriers['participant_phenotypic_sex'] != 'Indeterminate']
FD_noncarriers = FD_noncarriers.drop_duplicates(['participant_id'])
FD_noncarriers = FD_noncarriers.loc[FD_noncarriers['participant_phenotypic_sex'] != 'Indeterminate']
ARVC_noncarriers = ARVC_noncarriers.drop_duplicates(['participant_id'])
ARVC_noncarriers = ARVC_noncarriers.loc[ARVC_noncarriers['participant_phenotypic_sex'] != 'Indeterminate']
LQTS_BS_noncarriers = LQTS_BS_noncarriers.drop_duplicates(['participant_id'])
LQTS_BS_noncarriers = LQTS_BS_noncarriers.loc[LQTS_BS_noncarriers['participant_phenotypic_sex'] != 'Indeterminate']

# Recreating non-carrier dictionary
disease_lists = [HCM_noncarriers, DCM_noncarriers, FD_noncarriers, ARVC_noncarriers, LQTS_BS_noncarriers]

for i in range(len(disease_lists)):
    print(len(disease_lists[i]))
    print('')

# Writing data to files
HCM_noncarriers.to_csv("HCM_noncarriers_demographics.txt", index = False)
DCM_noncarriers.to_csv("DCM_noncarriers_demographics.txt", index = False)
FD_noncarriers.to_csv("FD_noncarriers_demographics.txt", index = False)
ARVC_noncarriers.to_csv("ARVC_noncarriers_demographics.txt", index = False)
LQTS_BS_noncarriers.to_csv("LQTS_BS_noncarriers_demographics.txt", index = False)
print("Data written to file!")

# Creating entire non-carrier dataframe
noncarriers = pd.concat(disease_lists)

# Creating non-carrier ID list
noncarrier_ids = noncarriers[['participant_id']]
noncarrier_ids = noncarrier_ids.drop_duplicates()
# Writing to file
noncarrier_ids.to_csv("non_carrier_ids.txt", index = False)
print("Data written to file!")

# # Building Demographics for Non-Carriers

# Dropping duplicates (if any)
HCM_noncarriers = HCM_noncarriers.drop_duplicates()
DCM_noncarriers = DCM_noncarriers.drop_duplicates()
FD_noncarriers = FD_noncarriers.drop_duplicates()
ARVC_noncarriers = ARVC_noncarriers.drop_duplicates()
LQTS_BS_noncarriers = LQTS_BS_noncarriers.drop_duplicates()

# Dropping duplicate rows per participant ID
HCM_noncarriers = HCM_noncarriers.drop_duplicates(['participant_id'])
DCM_noncarriers = DCM_noncarriers.drop_duplicates(['participant_id'])
FD_noncarriers = FD_noncarriers.drop_duplicates(['participant_id'])
ARVC_noncarriers = ARVC_noncarriers.drop_duplicates(['participant_id'])
LQTS_BS_noncarriers = LQTS_BS_noncarriers.drop_duplicates(['participant_id'])

# Dropping oarticipants where sex is indeterminate
HCM_noncarriers = HCM_noncarriers.loc[HCM_noncarriers['participant_phenotypic_sex'] != 'Indeterminate']
DCM_noncarriers = DCM_noncarriers.loc[DCM_noncarriers['participant_phenotypic_sex'] != 'Indeterminate']
FD_noncarriers = FD_noncarriers.loc[FD_noncarriers['participant_phenotypic_sex'] != 'Indeterminate']
ARVC_noncarriers = ARVC_noncarriers.loc[ARVC_noncarriers['participant_phenotypic_sex'] != 'Indeterminate']
LQTS_BS_noncarriers = LQTS_BS_noncarriers.loc[LQTS_BS_noncarriers['participant_phenotypic_sex'] != 'Indeterminate']

# Creating new dictionary with dataframes
disease_lists = {'HCM':HCM_noncarriers, 'DCM':DCM_noncarriers, 'FD':FD_noncarriers,
                 'ARVC':ARVC_noncarriers, 'LQTS & BS':LQTS_BS_noncarriers}

# Finding totals and any duplicates
for key in disease_lists:
    duplicates = disease_lists[key][disease_lists[key].duplicated(keep = False)]
    print(key)
    print(len(disease_lists[key]))
    print(len(duplicates))
    print('')

# Identifying gender counts
for key in disease_lists:
    print(key)
    print(disease_lists[key].participant_phenotypic_sex.value_counts())
    print('')

# Identifying ethnicity counts
for key in disease_lists:
    print(key)
    print(disease_lists[key].participant_ethnic_category_simple.value_counts().rename_axis('Ethnicity').reset_index(name='count'))
    print('')

# Identifying age category counts
for key in disease_lists:
    print(key)
    print(disease_lists[key].age_category.value_counts().rename_axis('Age Category').reset_index(name='count'))
    print('')

# Identifying ancestry counts
for key in disease_lists:
    print(key)
    print(disease_lists[key].ancestry.value_counts().rename_axis('Ancestry').reset_index(name='count'))
    print('')

# Identifying median height
for key in disease_lists:
    print(key, 'Median Height:')
    print(disease_lists[key]['height'].median())
    print('')

# Identifying median weight
for key in disease_lists:
    print(key, 'Median Weight:')
    print(disease_lists[key]['weight'].median())
    print('')

# Identifying age of onset category counts
for key in disease_lists:
    print(key)
    print(disease_lists[key].onset_category.value_counts().rename_axis('Age of Onset Category').reset_index(name='count'))
    print('')

# Creating non-carrier demographic table
# Total non-carriers
noncarrier_demographic_table = pd.DataFrame([[len(HCM_noncarriers.participant_id.unique()),
                                              len(DCM_noncarriers.participant_id.unique()),
                                              len(FD_noncarriers.participant_id.unique()),
                                            len(ARVC_noncarriers.participant_id.unique()), 
                                            len(LQTS_BS_noncarriers.participant_id.unique())]],
                                         index= ['Total Probands and Relatives'],
                                         columns = ['HCM', 'DCM', 'FD', 'ARVC', 'LQTS_BS'])

# Gender
gender_table = pd.DataFrame({'HCM': [41307, 36478],
                             'DCM': [41366, 36525],
                             'FD': [41427, 36594],
                             'ARVC': [41288, 36506],
                             'LQTS_BS': [41248, 36441]}, 
                            index=pd.MultiIndex.from_tuples(
                                [('Gender','Females'), ('Gender','Males')], 
                                names=['Demographic', 'Subgroup']))

noncarrier_demographic_table = pd.concat([noncarrier_demographic_table, gender_table])

# Ethnicity
ethnicity_table = pd.DataFrame({'HCM': [54535, 12185, 6428, 1847, 1437, 1353],
                                'DCM': [54637, 12203, 6436, 1832, 1427, 1356],
                                'FD': [54698, 12226, 6437, 1858, 1440, 1362],
                             'ARVC': [54530, 12195, 6430, 1847, 1432, 1360],
                             'LQTS_BS': [54479, 12176, 6407, 1841, 1434, 1352]}, 
                            index=pd.MultiIndex.from_tuples(
                                [('Ethnicity','White'), ('Ethnicity','Not known / Missing'), 
                                 ('Ethnicity','Asian or Asian British'), ('Ethnicity','Black or Black British'),
                                 ('Ethnicity','Mixed'), ('Ethnicity','Other')], 
                                names=['Demographic', 'Subgroup']))


noncarrier_demographic_table = pd.concat([noncarrier_demographic_table, ethnicity_table])

# Ancestry
ancestry_table = pd.DataFrame({'HCM': [61137, 16648],
                               'DCM': [61252, 16639],
                               'FD': [61327, 16694],
                             'ARVC': [61144, 16650],
                             'LQTS_BS': [61092, 16597]}, 
                            index=pd.MultiIndex.from_tuples(
                                [('Ancestry','European'), ('Ancestry','Other')], 
                                names=['Demographic', 'Subgroup']))

noncarrier_demographic_table = pd.concat([noncarrier_demographic_table, ancestry_table])

# Median height
height_table = pd.DataFrame({'HCM': [1.60],
                             'DCM': [1.60],
                             'FD': [1.60],
                             'ARVC': [1.60],
                             'LQTS_BS': [1.60]}, 
                            index=pd.MultiIndex.from_tuples(
                                [('Height','Median (m)')], 
                                names=['Demographic', 'Subgroup']))

noncarrier_demographic_table = pd.concat([noncarrier_demographic_table, height_table])

# Median weight
weight_table = pd.DataFrame({'HCM': [63.2],
                             'DCM': [63.3],
                             'FD': [63.2],
                             'ARVC': [63.2],
                             'LQTS_BS': [63.1]}, 
                            index=pd.MultiIndex.from_tuples(
                                [('Weight','Median (kg)')], 
                                names=['Demographic', 'Subgroup']))

noncarrier_demographic_table = pd.concat([noncarrier_demographic_table, weight_table])

# Age
age_table = pd.DataFrame({'HCM': [16332, 29013, 26072, 6284, '<5'],
                          'DCM': [16332, 29053, 26133, 6290, '<5'],
                          'FD': [16370, 29113, 26158, 6296, '<5'],
                             'ARVC': [16332, 29020, 26075, 6284, '<5'],
                             'LQTS_BS': [16301, 28981, 26053, 6271, '<5']}, 
                            index=pd.MultiIndex.from_tuples(
                                [('Age','0-25'), ('Age','25-50'), 
                                 ('Age','50-75'), ('Age','75-100'),
                                 ('Age','100-115')], 
                                names=['Demographic', 'Subgroup']))

noncarrier_demographic_table = pd.concat([noncarrier_demographic_table, age_table])

# Age of Onset
onset_table = pd.DataFrame({'HCM': [14630, 5010, 2591, 2461, 2474, 2405, 2636],
                            'DCM': [14637, 5015, 2600, 2474, 2494, 2415, 2645],
                            'FD': [14669, 5026, 2608, 2478, 2496, 2422, 2650],
                             'ARVC': [14627, 5012, 2599, 2466, 2488, 2411, 2647],
                             'LQTS_BS': [14618, 5001, 2592, 2465, 2488, 2407, 2643]}, 
                            index=pd.MultiIndex.from_tuples(
                                [('Age of Onset','Prenatal'), ('Age of Onset','0-10'), ('Age of Onset','10-20'), 
                                 ('Age of Onset','20-30'), ('Age of Onset','30-40'),
                                 ('Age of Onset','40-50'), ('Age of Onset', '50-115')], 
                                names=['Demographic', 'Subgroup']))

noncarrier_demographic_table = pd.concat([noncarrier_demographic_table, onset_table])

# Renaming columns
noncarrier_demographic_table = noncarrier_demographic_table.rename(columns = {'HCM':'HCM_non_carriers',
                                                                              'DCM':'DCM_non_carriers',
                                                                              'FD':'FD_non_carriers',
                                                                          'ARVC':'ARVC_non_carriers',
                                                                          'LQTS_BS':'LQTS_BS_non_carriers'})

# Adding column indicating label for each row
noncarrier_demographic_table['category'] = ['Total Probands and Relatives' ,'(Gender, Females)', '(Gender, Males)',
                                            '(Ethnicity, White)','(Ethnicity, Not known / Missing)',
                                            '(Ethnicity, Asian or Asian British)','(Ethnicity, Black or Black British)',
                                            '(Ethnicity, Mixed)','(Ethnicity, Other)','(Ancestry, European)',
                                            '(Ancestry, Other)','(Death, Count)','(Height, Median (m))',
                                            '(Weight, Median (kg))','(Age, 0-25)','(Age, 25-50)','(Age, 50-75)',
                                            '(Age, 75-100)','(Age, 100-115)','(Age of Onset, Prenatal)',
                                            '(Age of Onset, 0-10)','(Age of Onset, 10-20)' ,'(Age of Onset, 20-30)',
                                            '(Age of Onset, 30-40)','(Age of Onset, 40-50)','(Age of Onset, 50-115)' ]

# Writing to file
noncarrier_demographic_table.to_csv("final_non_carrier_demographics_table.txt", index = False)
print("Data written to file!")


# DEMOGRAPHICS STATISTICAL TESTS

# Loading in carrier and non-carrier IDs
carrier_ids = pd.read_csv("carrier_ids.txt", sep = ",")
non_carrier_ids = pd.read_csv("non_carrier_ids.txt", sep = ",")

# Creating all participant ID list
participant_ids = pd.concat([carrier_ids,non_carrier_ids])
participant_ids = participant_ids.drop_duplicates()

# All participant IDs file
participant_ids.to_csv("all_participant_ids.txt", index = False)
print("Data written to file!")

# Loading in carrier and non-carrier demographic tables
carriers = pd.read_csv("final_carrier_demographics_table.txt", sep = ",")
non_carriers = pd.read_csv("final_non_carrier_demographics_table.txt", sep = ",")

# Merging carrier and non-carrier tables
demographics = pd.merge(carriers, non_carriers, on='category', how='left')

# Reordering columns
demographics = demographics[['category', 'HCM_non_carriers', 'HCM_carriers',
                             'DCM_non_carriers', 'DCM_carriers',
                             'FD_non_carriers', 'FD_carriers',
                            'ARVC_non_carriers', 'ARVC_carriers', 
                            'LQTS_BS_non_carriers', 'LQTS_BS_carriers']]

# Writing to file
demographics.to_csv("final_carriers_vs_non_carriers_numbers.txt", index = False)
print("Data written to file!")

# # Statistics
# Loading in carrier and non-carrier files
HCM_carriers = pd.read_csv("HCM_carrier_demographics.txt", sep = ",")
DCM_carriers = pd.read_csv("DCM_carrier_demographics.txt", sep = ",")
FD_carriers = pd.read_csv("FD_carrier_demographics.txt", sep = ",")
ARVC_carriers = pd.read_csv("ARVC_carrier_demographics.txt", sep = ",")
LQTS_BS_carriers = pd.read_csv("LQTS_BS_carrier_demographics.txt", sep = ",")
HCM_non_carriers = pd.read_csv("HCM_noncarriers_demographics.txt", sep = ",")
DCM_non_carriers = pd.read_csv("DCM_noncarriers_demographics.txt", sep = ",")
FD_non_carriers = pd.read_csv("FD_noncarriers_demographics.txt", sep = ",")
ARVC_non_carriers = pd.read_csv("ARVC_noncarriers_demographics.txt", sep = ",")
LQTS_BS_non_carriers = pd.read_csv("LQTS_BS_noncarriers_demographics.txt", sep = ",")

# Merging carrier and non-carrier files for each disease, adding column indicating carrier status
HCM = pd.merge(HCM_carriers, HCM_non_carriers, how='outer')
HCM['carrier_status'] = np.where(HCM['participant_id'].isin(HCM_carriers.participant_id), 1, 0) # carrier = 1
DCM = pd.merge(DCM_carriers, DCM_non_carriers, how='outer')
DCM['carrier_status'] = np.where(DCM['participant_id'].isin(DCM_carriers.participant_id), 1, 0) 
FD = pd.merge(FD_carriers, FD_non_carriers, how='outer')
FD['carrier_status'] = np.where(FD['participant_id'].isin(FD_carriers.participant_id), 1, 0) 
ARVC = pd.merge(ARVC_carriers, ARVC_non_carriers, how='outer')
ARVC['carrier_status'] = np.where(ARVC['participant_id'].isin(ARVC_carriers.participant_id), 1, 0)
LQTS_BS = pd.merge(LQTS_BS_carriers, LQTS_BS_non_carriers, how='outer')
LQTS_BS['carrier_status'] = np.where(LQTS_BS['participant_id'].isin(LQTS_BS_carriers.participant_id), 1, 0)

# Setting not died status to 0 instead of null
HCM["died"] = HCM["died"].fillna(0)
DCM["died"] = DCM["died"].fillna(0)
FD["died"] = FD["died"].fillna(0)
ARVC["died"] = ARVC["died"].fillna(0)
LQTS_BS["died"] = LQTS_BS["died"].fillna(0)

# Ensuring all dataframes have one row per participant and participants with missing sex removed
HCM = HCM.drop_duplicates(['participant_id', 'gene'])
HCM = HCM.loc[HCM['participant_phenotypic_sex'] != 'Indeterminate']
DCM = DCM.drop_duplicates(['participant_id', 'gene'])
DCM = DCM.loc[DCM['participant_phenotypic_sex'] != 'Indeterminate']
FD = FD.drop_duplicates(['participant_id', 'gene'])
FD = FD.loc[FD['participant_phenotypic_sex'] != 'Indeterminate']
ARVC = ARVC.drop_duplicates(['participant_id', 'gene'])
ARVC = ARVC.loc[ARVC['participant_phenotypic_sex'] != 'Indeterminate']
LQTS_BS = LQTS_BS.drop_duplicates(['participant_id', 'gene'])
LQTS_BS = LQTS_BS.loc[LQTS_BS['participant_phenotypic_sex'] != 'Indeterminate']

import scipy.stats as stats

# # Phenotypic Sex

# Chi-Squared test for HCM gender between carriers and non-carriers
HCM_crosstab_sex = pd.crosstab(HCM['participant_phenotypic_sex'], HCM['carrier_status'])
stats.chi2_contingency(HCM_crosstab_sex)

# Chi-Squared test for DCM gender between carriers and non-carriers
DCM_crosstab_sex = pd.crosstab(DCM['participant_phenotypic_sex'], DCM['carrier_status'])
stats.chi2_contingency(DCM_crosstab_sex)

# Chi-Squared test for FD gender between carriers and non-carriers
FD_crosstab_sex = pd.crosstab(FD['participant_phenotypic_sex'], FD['carrier_status'])
stats.chi2_contingency(FD_crosstab_sex)

# Chi-Squared test for ARVC gender between carriers and non-carriers
ARVC_crosstab_sex = pd.crosstab(ARVC['participant_phenotypic_sex'], ARVC['carrier_status'])
stats.chi2_contingency(ARVC_crosstab_sex)

# Chi-Squared test for LQTS gender between carriers and non-carriers
LQTS_BS_crosstab_sex = pd.crosstab(LQTS_BS['participant_phenotypic_sex'], LQTS_BS['carrier_status'])
stats.chi2_contingency(LQTS_BS_crosstab_sex)

# # Age

# Chi-Squared test for HCM age between carriers and non-carriers
HCM_crosstab_age = pd.crosstab(HCM['age_category'], HCM['carrier_status'])
stats.chi2_contingency(HCM_crosstab_age)

# Chi-Squared test for DCM age between carriers and non-carriers
DCM_crosstab_age = pd.crosstab(DCM['age_category'], DCM['carrier_status'])
stats.chi2_contingency(DCM_crosstab_age)

# Chi-Squared test for FD age between carriers and non-carriers
FD_crosstab_age = pd.crosstab(FD['age_category'], FD['carrier_status'])
stats.chi2_contingency(FD_crosstab_age)

# Chi-Squared test for ARVC age between carriers and non-carriers
ARVC_crosstab_age = pd.crosstab(ARVC['age_category'], ARVC['carrier_status'])
stats.chi2_contingency(ARVC_crosstab_age)

# Chi-Squared test for LQTS age between carriers and non-carriers
LQTS_BS_crosstab_age = pd.crosstab(LQTS_BS['age_category'], LQTS_BS['carrier_status'])
stats.chi2_contingency(LQTS_BS_crosstab_age)

# # Ethnic Category

# Chi-Squared test for HCM ethnicity between carriers and non-carriers
HCM_crosstab_ethnic = pd.crosstab(HCM['participant_ethnic_category_simple'], HCM['carrier_status'])
stats.chi2_contingency(HCM_crosstab_ethnic)

# Chi-Squared test for DCM ethnicity between carriers and non-carriers
DCM_crosstab_ethnic = pd.crosstab(DCM['participant_ethnic_category_simple'], DCM['carrier_status'])
stats.chi2_contingency(DCM_crosstab_ethnic)

# Chi-Squared test for FD ethnicity between carriers and non-carriers
FD_crosstab_ethnic = pd.crosstab(FD['participant_ethnic_category_simple'], FD['carrier_status'])
stats.chi2_contingency(FD_crosstab_ethnic)

# Chi-Squared test for ARVC ethnicity between carriers and non-carriers
ARVC_crosstab_ethnic = pd.crosstab(ARVC['participant_ethnic_category_simple'], ARVC['carrier_status'])
stats.chi2_contingency(ARVC_crosstab_ethnic)

# Chi-Squared test for LQTS ethnicity between carriers and non-carriers
LQTS_BS_crosstab_ethnic = pd.crosstab(LQTS_BS['participant_ethnic_category_simple'], LQTS_BS['carrier_status'])
stats.chi2_contingency(LQTS_BS_crosstab_ethnic)

# # Ancestry

# Chi-Squared test for HCM ancestry between carriers and non-carriers
HCM_crosstab_ancestry = pd.crosstab(HCM['ancestry'], HCM['carrier_status'])
stats.chi2_contingency(HCM_crosstab_ancestry)

# Chi-Squared test for DCM ancestry between carriers and non-carriers
DCM_crosstab_ancestry = pd.crosstab(DCM['ancestry'], DCM['carrier_status'])
stats.chi2_contingency(DCM_crosstab_ancestry)

# Chi-Squared test for FD ancestry between carriers and non-carriers
FD_crosstab_ancestry = pd.crosstab(FD['ancestry'], FD['carrier_status'])
stats.chi2_contingency(FD_crosstab_ancestry)

# Chi-Squared test for ARVC ancestry between carriers and non-carriers
ARVC_crosstab_ancestry = pd.crosstab(ARVC['ancestry'], ARVC['carrier_status'])
stats.chi2_contingency(ARVC_crosstab_ancestry)

# Chi-Squared test for LQTS ancestry between carriers and non-carriers
LQTS_BS_crosstab_ancestry = pd.crosstab(LQTS_BS['ancestry'], LQTS_BS['carrier_status'])
stats.chi2_contingency(LQTS_BS_crosstab_ancestry)


# # Age of Onset

# Chi-Squared test for HCM onset between carriers and non-carriers
HCM_crosstab_onset = pd.crosstab(HCM['onset_category'], HCM['carrier_status'])
stats.chi2_contingency(HCM_crosstab_onset)

# Chi-Squared test for DCM onset between carriers and non-carriers
DCM_crosstab_onset = pd.crosstab(DCM['onset_category'], DCM['carrier_status'])
stats.chi2_contingency(DCM_crosstab_onset)

# Chi-Squared test for FD onset between carriers and non-carriers
FD_crosstab_onset = pd.crosstab(FD['onset_category'], FD['carrier_status'])
stats.chi2_contingency(FD_crosstab_onset)

# Chi-Squared test for ARVC onset between carriers and non-carriers
ARVC_crosstab_onset = pd.crosstab(ARVC['onset_category'], ARVC['carrier_status'])
stats.chi2_contingency(ARVC_crosstab_onset)

# Chi-Squared test for LQTS onset between carriers and non-carriers
LQTS_BS_crosstab_onset = pd.crosstab(LQTS_BS['onset_category'], LQTS_BS['carrier_status'])
stats.chi2_contingency(LQTS_BS_crosstab_onset)

# ## Creating t-test function
# Source:
# https://machinelearningmastery.com/how-to-code-the-students-t-test-from-scratch-in-python/

from math import sqrt
from numpy import mean
from scipy.stats import sem
from scipy.stats import t

def independent_ttest(data1, data2, alpha):
    # calculate means
    mean1, mean2 = mean(data1), mean(data2)
    # calculate standard errors
    se1, se2 = sem(data1), sem(data2)
    # standard error on the difference between the samples
    sed = sqrt(se1**2.0 + se2**2.0)
    # calculate the t statistic
    t_stat = (mean1 - mean2) / sed
    # degrees of freedom
    df = len(data1) + len(data2) - 2
    # calculate the critical value
    cv = t.ppf(1.0 - alpha, df)
    # calculate the p-value
    p = (1.0 - t.cdf(abs(t_stat), df)) * 2.0
    # return everything
    return t_stat, df, cv, p

# # Height

# Applying t-test to assess difference between mean heights of HCM carriers and non-carriers
# calculate the t test
alpha = 0.05
t_stat, df, cv, p = independent_ttest(HCM[HCM['carrier_status'] == 1]['height'], 
                                      HCM[HCM['carrier_status'] == 0]['height'], alpha)
print('t=%.3f, df=%d, cv=%.3f, p=%.3f' % (t_stat, df, cv, p))

# interpret via p-value
if p > alpha:
    print('Accept null hypothesis that the means are equal.')
else:
    print('Reject the null hypothesis that the means are equal.')

# Applying t-test to assess difference between mean heights of DCM carriers and non-carriers
t_stat, df, cv, p = independent_ttest(DCM[DCM['carrier_status'] == 1]['height'], 
                                      DCM[DCM['carrier_status'] == 0]['height'], alpha)
print('t=%.3f, df=%d, cv=%.3f, p=%.3f' % (t_stat, df, cv, p))
# interpret via p-value
if p > alpha:
    print('Accept null hypothesis that the means are equal.')
else:
    print('Reject the null hypothesis that the means are equal.')

# Applying t-test to assess difference between mean heights of FD carriers and non-carriers
t_stat, df, cv, p = independent_ttest(FD[FD['carrier_status'] == 1]['height'], 
                                      FD[FD['carrier_status'] == 0]['height'], alpha)
print('t=%.3f, df=%d, cv=%.3f, p=%.3f' % (t_stat, df, cv, p))
# interpret via p-value
if p > alpha:
    print('Accept null hypothesis that the means are equal.')
else:
    print('Reject the null hypothesis that the means are equal.')

# Applying t-test to assess difference between mean heights of ARVC carriers and non-carriers
t_stat, df, cv, p = independent_ttest(ARVC[ARVC['carrier_status'] == 1]['height'], 
                                      ARVC[ARVC['carrier_status'] == 0]['height'], alpha)
print('t=%.3f, df=%d, cv=%.3f, p=%.3f' % (t_stat, df, cv, p))
# interpret via p-value
if p > alpha:
    print('Accept null hypothesis that the means are equal.')
else:
    print('Reject the null hypothesis that the means are equal.')

# Applying t-test to assess difference between mean heights of LQTS carriers and non-carriers
t_stat, df, cv, p = independent_ttest(LQTS_BS[LQTS_BS['carrier_status'] == 1]['height'], 
                                      LQTS_BS[LQTS_BS['carrier_status'] == 0]['height'], alpha)
print('t=%.3f, df=%d, cv=%.3f, p=%.3f' % (t_stat, df, cv, p))
# interpret via p-value
if p > alpha:
    print('Accept null hypothesis that the means are equal.')
else:
    print('Reject the null hypothesis that the means are equal.')

# # Weight

# Applying t-test to assess difference between mean weights of HCM carriers and non-carriers
t_stat, df, cv, p = independent_ttest(HCM[HCM['carrier_status'] == 1]['weight'], 
                                      HCM[HCM['carrier_status'] == 0]['weight'], alpha)
print('t=%.3f, df=%d, cv=%.3f, p=%.3f' % (t_stat, df, cv, p))
# interpret via p-value
if p > alpha:
    print('Accept null hypothesis that the means are equal.')
else:
    print('Reject the null hypothesis that the means are equal.')

# Applying t-test to assess difference between mean weights of DCM carriers and non-carriers
t_stat, df, cv, p = independent_ttest(DCM[DCM['carrier_status'] == 1]['weight'], 
                                      DCM[DCM['carrier_status'] == 0]['weight'], alpha)
print('t=%.3f, df=%d, cv=%.3f, p=%.3f' % (t_stat, df, cv, p))
# interpret via p-value
if p > alpha:
    print('Accept null hypothesis that the means are equal.')
else:
    print('Reject the null hypothesis that the means are equal.')

# Applying t-test to assess difference between mean weights of FD carriers and non-carriers
t_stat, df, cv, p = independent_ttest(FD[FD['carrier_status'] == 1]['weight'], 
                                      FD[FD['carrier_status'] == 0]['weight'], alpha)
print('t=%.3f, df=%d, cv=%.3f, p=%.3f' % (t_stat, df, cv, p))
# interpret via p-value
if p > alpha:
    print('Accept null hypothesis that the means are equal.')
else:
    print('Reject the null hypothesis that the means are equal.')

# Applying t-test to assess difference between mean weights of ARVC carriers and non-carriers
t_stat, df, cv, p = independent_ttest(ARVC[ARVC['carrier_status'] == 1]['weight'], 
                                      ARVC[ARVC['carrier_status'] == 0]['weight'], alpha)
print('t=%.3f, df=%d, cv=%.3f, p=%.3f' % (t_stat, df, cv, p))
# interpret via p-value
if p > alpha:
    print('Accept null hypothesis that the means are equal.')
else:
    print('Reject the null hypothesis that the means are equal.')

# Applying t-test to assess difference between mean weights of LQTS carriers and non-carriers
t_stat, df, cv, p = independent_ttest(LQTS_BS[LQTS_BS['carrier_status'] == 1]['weight'], 
                                      LQTS_BS[LQTS_BS['carrier_status'] == 0]['weight'], alpha)
print('t=%.3f, df=%d, cv=%.3f, p=%.3f' % (t_stat, df, cv, p))
# interpret via p-value
if p > alpha:
    print('Accept null hypothesis that the means are equal.')
else:
    print('Reject the null hypothesis that the means are equal.')