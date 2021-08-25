import labkey
import pandas as pd
import numpy as np

# PENETRANCE
# Loading in table rare_diseases_participant_disease (info on recruitment participant disease)
labkey_server = "labkey-embassy.gel.zone" 
project_name = "main-programme/main-programme_v11_2020-12-17"  
context_path = "labkey"  
schema_name = "lists"  
query_name = "rare_diseases_participant_disease"  
server_context = labkey.utils.create_server_context(
    labkey_server, project_name, context_path, use_ssl=True
)
disease_results = labkey.query.select_rows(server_context, schema_name, query_name, max_rows=200000)
disease_results = pd.DataFrame(disease_results["rows"])

# Subsetting for column with recruitment disease
recruitment_diseases = disease_results[['participant_id','specific_disease']]

# loading in data on carriers of variants
carriers = pd.read_csv("full_carriers_info.txt", sep = ",")
carriers.head()

carriers = carriers.drop_duplicates(['participant_id', 'gene'])
carriers = carriers.loc[carriers['participant_phenotypic_sex'] != 'Indeterminate']

# dropping disease from this file to re-do and make sure correct
carriers = carriers.drop(columns=['specific_disease'])

# Merging carrier information with the disease they were recruited for
data = pd.merge(carriers, recruitment_diseases, how='left', on='participant_id')

# Renaming column for recruitment disease
data = data.rename(columns = {'specific_disease':'recruitment_disease'})

# Subsetting for wanted information
penetrance = data[['participant_id', 
                   'participant_type', 
                   'gene', 
                   'ID', 
                   'genotype', 
                   'source', 
                   'recruitment_disease']]

# create a list of conditions (genes for each disease)
conditions = [
    (penetrance['gene'] == 'MYBPC3') | 
    (penetrance['gene'] == 'MYH7') |
    (penetrance['gene'] == 'TNNI3') |
    (penetrance['gene'] == 'TPM1') |
    (penetrance['gene'] == 'MYL3') |
    (penetrance['gene'] == 'ACTC1') |
    (penetrance['gene'] == 'MYL2') |
    (penetrance['gene'] == 'PRKAG2'),
    (penetrance['gene'] == 'TNNT2'),
    (penetrance['gene'] == 'GLA'),
    (penetrance['gene'] == 'LMNA'),
    (penetrance['gene'] == 'PKP2') | (penetrance['gene'] == 'DSP') | (penetrance['gene'] == 'DSC2') | (penetrance['gene'] == 'TMEM43') | (penetrance['gene'] == 'DSG2'),
    (penetrance['gene'] == 'KCNQ1') | (penetrance['gene'] == 'KCNH2') | (penetrance['gene'] == 'SCN5A')    
     ]

# create a list of the diseases to assign for each condition
values = ['Hypertrophic Cardiomyopathy', 
          'Dilated Cardiomyopathy / Left Ventricular Noncompaction', 
          'Fabry Disease (Hypertrophic Cardiomyopathy)', 
          'Dilated Cardiomyopathy',
          'Arrhythmogenic Right Ventricular Cardiomyopathy', 
          'Long QT / Brugada Syndrome']

# create a new column and use np.select to assign values to it using our lists as arguments
penetrance['variant_specific_disease'] = np.select(conditions, values)

# Looking at possible recruitment diseases
penetrance['recruitment_disease'].unique()

# ### My terms
# 'Hypertrophic Cardiomyopathy', 
# 'Long QT syndrome',
# 'Brugada syndrome', 
# 'Arrhythmogenic Right Ventricular Cardiomyopathy',
# 'Dilated Cardiomyopathy',
# 'Left Ventricular Noncompaction Cardiomyopathy', 
# 'Dilated Cardiomyopathy and conduction defects',
# 'Dilated Cardiomyopathy (DCM)',
# 'Mucopolysaccharideosis, Gaucher, Fabry',
# 'Mucopolysaccharideosis -  Gaucher -  Fabry',
# 'hypertrophic cardiomyopathy' 
# ###

# Creating my recruited disease column with my terms for diseases to later compare to variant disease
penetrance['my_recruitment_disease'] = np.nan

# Adding in correct terms and replacing with standard term
penetrance['my_recruitment_disease'] = np.where((penetrance['recruitment_disease'] == 'Hypertrophic Cardiomyopathy') | 
                                                (penetrance['recruitment_disease'] == 'hypertrophic cardiomyopathy') , 
                                                'Hypertrophic Cardiomyopathy', 
                                                penetrance['my_recruitment_disease']) 
penetrance['my_recruitment_disease'] = np.where(penetrance['recruitment_disease'] == 'Left Ventricular Noncompaction Cardiomyopathy', 
                                                'Dilated Cardiomyopathy / Left Ventricular Noncompaction', 
                                                penetrance['my_recruitment_disease']) 
penetrance['my_recruitment_disease'] = np.where((penetrance['recruitment_disease'] == 'Mucopolysaccharideosis, Gaucher, Fabry') |
                                                (penetrance['recruitment_disease'] == 'Mucopolysaccharideosis -  Gaucher -  Fabry'), 
                                                'Fabry Disease (Hypertrophic Cardiomyopathy)', 
                                                penetrance['my_recruitment_disease']) 
penetrance['my_recruitment_disease'] = np.where((penetrance['recruitment_disease'] == 'Dilated Cardiomyopathy') |
                                                (penetrance['recruitment_disease'] == 'Dilated Cardiomyopathy and conduction defects') |
                                                (penetrance['recruitment_disease'] == 'Dilated Cardiomyopathy (DCM)'), 
                                                'Dilated Cardiomyopathy', 
                                                penetrance['my_recruitment_disease']) 
penetrance['my_recruitment_disease'] = np.where(penetrance['recruitment_disease'] == 'Arrhythmogenic Right Ventricular Cardiomyopathy', 
                                                'Arrhythmogenic Right Ventricular Cardiomyopathy', 
                                                penetrance['my_recruitment_disease']) 
penetrance['my_recruitment_disease'] = np.where((penetrance['recruitment_disease'] == 'Long QT syndrome') |
                                                (penetrance['recruitment_disease'] == 'Brugada syndrome'), 
                                                'Long QT / Brugada Syndrome', 
                                                penetrance['my_recruitment_disease']) 

# Adding in recruited disease to my recruitement disease column
penetrance['my_recruitment_disease'] = np.where((penetrance['my_recruitment_disease'] == 'nan'), 
                                                penetrance['recruitment_disease'], 
                                                penetrance['my_recruitment_disease']) 

# Adding column indicating if variant specific disease == recruited disease
# 0 = no, 1 = yes
penetrance['recruited_for_variant_disease'] = np.nan
penetrance['recruited_for_variant_disease'] = np.where((penetrance['my_recruitment_disease'] == penetrance['variant_specific_disease']), 
                                                       1,
                                                       penetrance['recruited_for_variant_disease'])

# If recruited for variant is null, input 0 
penetrance['recruited_for_variant_disease'] = np.where((penetrance['recruited_for_variant_disease'].isnull()), 
                                                       0, 
                                                       penetrance['recruited_for_variant_disease']) 


# Observing if any duplicates 
duplicates = penetrance[penetrance.duplicated(['participant_id', 'gene', 'ID', 
                                               'variant_specific_disease', 'my_recruitment_disease'], keep = False)]

# Dropping duplicates
penetrance = penetrance.drop_duplicates(['participant_id', 'gene', 'ID', 
                                               'variant_specific_disease', 'my_recruitment_disease'])

# Subsetting for wanted information
penetrance = penetrance[['participant_id', 
                   'participant_type', 
                   'gene', 
                   'ID', 
                   'genotype', 
                   'variant_specific_disease', 
                   'my_recruitment_disease',
                   'recruited_for_variant_disease']]

# Writing to file
penetrance.to_csv("penetrance.txt", index = False)
print("Data written to file!")

# loading in data 
carriers = pd.read_csv("full_carriers_info.txt", sep = ",")
penetrance = pd.read_csv("penetrance.txt", sep = ",")
icd10 = pd.read_csv("icd10.txt.gz", sep = "\t")
participants = pd.read_csv("participants_in_gvcf.txt", sep = "\t")

# List of genes of interest
genes = ['MYBPC3','MYH7','TNNT2','TNNI3','TPM1','MYL3','ACTC1','PRKAG2','GLA','MYL2',
         'LMNA','PKP2','DSP','DSC2','TMEM43','DSG2','KCNQ1','KCNH2','SCN5A']

# Making sure dataframe only contains genes of interest
penetrance = penetrance.loc[penetrance['gene'].isin(genes)]

# Participants with more than one row
duplicates = penetrance[penetrance.duplicated(['participant_id', 'gene'], keep = False)]
duplicates.head(20) # participants with more than 1 recruitment disease

# dropping duplicates with 2 recruitment diseases 
penetrance = penetrance.drop_duplicates(['participant_id', 'gene'])

# Participants with more than one row (carrying more than one variant)
duplicates = penetrance[penetrance.duplicated(['participant_id'], keep = False)]
print(len(duplicates.participant_id.unique()))

# Merging ICD-10 information to carriers
total_penetrance = pd.merge(penetrance, icd10, on='participant_id', how='left')

# Subsetting columns of interest
all_penetrance = total_penetrance[['participant_id', 'participant_type', 'gene', 
                                     'variant_specific_disease', 'my_recruitment_disease',
                                     'recruited_for_variant_disease', 'icd10', 'date', 'table', 
                                     'column_name', 'QC_errors', 'QC_flags' ]]

# Reasons for errors in ICD-10 information
all_penetrance.QC_flags.unique()

# Where QC_errors, setting ICD-10 code to null value as dont want to include in analysis but want to keep
# carrier in numbers
all_penetrance['icd10'] = np.where((all_penetrance.QC_errors != 0), 
                                   'ICD-NA', 
                                   all_penetrance['icd10'])

# Tidying dataframe
all_penetrance = all_penetrance.reset_index()
all_penetrance = all_penetrance.drop(columns='index')

# Tables where ICD-10 info came from
all_penetrance.table.unique()

# Dropping columns that are no longer needed
all_penetrance = all_penetrance.drop(columns=['date', 'table', 'column_name', 'QC_errors', 'QC_flags'])

# Splitting out into diseases
HCM = all_penetrance.loc[(all_penetrance['gene'] == 'MYBPC3') | 
    (all_penetrance['gene'] == 'MYH7') |
    (all_penetrance['gene'] == 'TNNI3') |
    (all_penetrance['gene'] == 'TPM1') |
    (all_penetrance['gene'] == 'MYL3') |
    (all_penetrance['gene'] == 'ACTC1') |
    (all_penetrance['gene'] == 'MYL2') |
    (all_penetrance['gene'] == 'PRKAG2')]
DCM = all_penetrance.loc[(all_penetrance['gene'] == 'TNNT2') | 
                            (all_penetrance['gene'] == 'LMNA')]
FD = all_penetrance.loc[all_penetrance['gene'] == 'GLA']
ARVC = all_penetrance.loc[(all_penetrance['gene'] == 'PKP2') |
                             (all_penetrance['gene'] == 'DSP') |
                             (all_penetrance['gene'] == 'DSC2') |
                             (all_penetrance['gene'] == 'TMEM43') |
                             (all_penetrance['gene'] == 'DSG2')]
LQTS_BS = all_penetrance.loc[(all_penetrance['gene'] == 'KCNQ1') |
                                (all_penetrance['gene'] == 'KCNH2') |
                                (all_penetrance['gene'] == 'SCN5A')]

HCM = HCM.sort_values('participant_id')
DCM = DCM.sort_values('participant_id')
FD = FD.sort_values('participant_id')
ARVC = ARVC.sort_values('participant_id')
LQTS_BS = LQTS_BS.sort_values('participant_id')

# Creating dictionary with disease dataframes
phenotype_dict = {'HCM': HCM,
                 'DCM': DCM,
                 'FD': FD,
                 'ARVC': ARVC,
                 'LQTS_BS':LQTS_BS}

for key in phenotype_dict:
    print(key, ':')
    print('Total carriers:', len(phenotype_dict[key]['participant_id'].unique()))


# ## Related ICD-10 codes for each pathogenic variant
MYBPC3 = ['I630', 'I631', 'I632', 'I633', 'I634', 'I635', 'I638', 'I639', 'I693', 'R55', 'I451', 'I447', 'G45', 
          'G458', 'G459', 'I460', 'I461', 'I462', 'I463', 'I464', 'I465', 'I466', 'I467', 'I468', 'I469', 'I313', 
          'I490', 'I517', 'R060', 'J81', 'R160', 'R18', 'G729', 'R073', 'R074']

MYH7 = ['I480', 'I481', 'I482', 'I483', 'I484', 'I489', 'I110', 'I130', 'I132', 'I50', 'R96', 'R960', 'R961']

TNNT2 = ['I48', 'I480', 'I481', 'I482', 'I483', 'I484', 'I489', 'I110', 'I130', 'I132', 'I50', 'I420']

TNNI3 = ['I48', 'I480', 'I481', 'I482', 'I483', 'I484', 'I489', 'I42']

TPM1 = ['I421', 'I422']

MYL3 = ['I421', 'I422', 'I460', 'I461', 'I462', 'I463', 'I464', 'I465', 'I466', 'I467', 'I468', 'I469', 'I490', 
        'I421', 'I422', 'I110', 'I130', 'I132', 'I50', 'I423', 'I425', 'R060', 'I460', 'I461', 'I469']

ACTC1 = ['I480', 'I481', 'I482', 'I483', 'I484', 'I489', 'I421', 'I422']

PRKAG2 = ['I447', 'I456', 'I48', 'I480', 'I481', 'I482', 'I483', 'I484', 'I489', 'R001', 'I44', 'I440', 'I441', 
          'I442', 'I443', 'I45', 'I421', 'I422', 'I456']

GLA = ['E752', 'H185', 'G45', 'I201', 'I208', 'I209', 'I480', 'I481', 'I482', 'I483', 'I484', 'I489', 'I200', 
       'I201', 'I208', 'I209', 'I252', 'I21', 'I22', 'I23', 'I241',  'I110', 'I130', 'I132', 'I50', 'R11', 'K529', 
       'A099', 'P783', 'E300', 'N18', 'N19', 'D60', 'D61', 'D62', 'D63', 'D64', 'D50', 'D51', 'D52', 'D53', 'D55', 
       'D56', 'D57', 'D58', 'D59', 'L744', 'R568', 'P90', 'R202', 'G571', 'R253', 'R101', 'R103', 'R104', 'R80', 
       'N391', 'N392', 'O121', 'O122', 'Q820', 'I890']

MYL2 = ['R42', 'I470', 'I472', 'I471', 'R002', 'I490', 'I421', 'I422', 'R060', 'R073', 'R074', 'I460', 
        'I461', 'I469']

LMNA = ['I420', 'I471', 'I48', 'I480', 'I481', 'I482', 'I483', 'I484', 'I489', 'I460', 'I461', 'I462', 'I463', 
        'I464', 'I465', 'I466', 'I467', 'I468', 'I469', 'I44', 'I440', 'I441', 'I442', 'I443', 'I45', 'I420', 
        'I426', 'I427', 'O903', 'I423', 'I424']

PKP2 = ['I428', 'R55', 'R002', 'I470', 'I472', 'I517', 'I460', 'I461', 'I469']

DSP = ['I428', 'I470', 'I472', 'I490', 'I110', 'I130', 'I132', 'I50', 'I493', 'I460', 'I461', 'I469']

DSC2 = ['I428', 'R55', 'R002', 'I470', 'I472', 'R060', 'Q841', 'I460', 'I461', 'I469']

TMEM43 = ['I428', 'I470', 'I472', 'R002', 'R42', 'I110', 'I130', 'I132', 'I50', 'I493', 'I447', 'I451', 'G241', 
          'G242', 'G248', 'G249', 'R073', 'R074', 'I460', 'I461', 'I469']

DSG2 = ['I428', 'I470', 'I472', 'R002', 'I493', 'I460', 'I461', 'I469']

KCNQ1 = ['I458', 'R55', 'I472', 'I490', 'I460', 'I461', 'I469']

KCNH2 = ['I458', 'R55', 'I472', 'I490', 'I460', 'I461', 'I469']

SCN5A = ['I458', 'R55', 'I472', 'I490', 'I460', 'I461', 'I469']

# Creating dictionary with related ICD-10 codes for each pathogenic variant
genes = {'MYBPC3': MYBPC3,
         'MYH7': MYH7,
         'TNNT2': TNNT2,
         'TNNI3': TNNI3,
         'TPM1': TPM1, 
         'MYL3': MYL3, 
         'ACTC1': ACTC1,
         'PRKAG2': PRKAG2,
         'GLA': GLA,
         'MYL2': MYL2,
         'LMNA': LMNA,
         'PKP2': PKP2,
         'DSP': DSP,
         'DSC2': DSC2,
         'TMEM43': TMEM43,
         'DSG2': DSG2,
         'KCNQ1': KCNQ1,
         'KCNH2': KCNH2,
         'SCN5A': SCN5A}

# Adding columns to phenotype DFs that indicates if ICD-10 is related to variant specific disease
for key in phenotype_dict:
    phenotype_dict[key]['variant_related_icd10'] = np.nan
    for k in genes:
        phenotype_dict[key]['variant_related_icd10'] = np.where((phenotype_dict[key]['gene'] == k) & 
                                                           (phenotype_dict[key]['icd10'].isin(genes[k])), 
                                                           1, 
                                                           phenotype_dict[key]['variant_related_icd10']) 
    phenotype_dict[key]['variant_related_icd10'] = np.where((phenotype_dict[key]['variant_related_icd10'].isnull()), 
                                                       0, 
                                                       phenotype_dict[key]['variant_related_icd10']) 

for key in phenotype_dict:
    phenotype_dict[key] = phenotype_dict[key].sort_values(['participant_id', 'icd10'])
    phenotype_dict[key] = phenotype_dict[key].reset_index()
    phenotype_dict[key] = phenotype_dict[key].drop(columns='index')

# Merging all separate phenotype DFs into one large dataframe, with one row per participant per ICD-10 code
dfs = [HCM,DCM,FD,ARVC,LQTS_BS]
df = pd.concat(dfs)
single_icd = pd.DataFrame({'code_count' : df.groupby(['participant_id', 
                                                      'gene',
                                                      'icd10']).size()}).reset_index()

# Cleaning dataframe
single_icd = single_icd.sort_values('participant_id')
single_icd = single_icd.reset_index()
single_icd = single_icd.drop(columns='index')

# Participants with duplicate rows 
duplicates = single_icd[single_icd.duplicated(keep = False)]
duplicates.head(5)

# Ensuring only genes of interest in table
info = single_icd.loc[(single_icd['gene'] == 'MYBPC3') | 
                                     (single_icd['gene'] == 'MYH7') |
                                     (single_icd['gene'] == 'TNNT2') |
                                     (single_icd['gene'] == 'TNNI3') |
                                     (single_icd['gene'] == 'TPM1') |
                                     (single_icd['gene'] == 'MYL3') |
                                     (single_icd['gene'] == 'ACTC1') |
                                     (single_icd['gene'] == 'PRKAG2') |
                                     (single_icd['gene'] == 'GLA') |
                                     (single_icd['gene'] == 'MYL2') |
                                     (single_icd['gene'] == 'LMNA') |
                           (single_icd['gene'] == 'PKP2') |
                           (single_icd['gene'] == 'DSP') |
                           (single_icd['gene'] == 'DSC2') |
                           (single_icd['gene'] == 'TMEM43') |
                           (single_icd['gene'] == 'DSG2') |
                           (single_icd['gene'] == 'KCNQ1') |
                           (single_icd['gene'] == 'KCNH2') |
                           (single_icd['gene'] == 'SCN5A')]

# Getting information back on recruited for disease, variant related ICD-10
info = pd.merge(info, df, on=['participant_id',
                                      'gene',
                                      'icd10'],
               how = 'left')

# Adding age and age of onset information
age_info = carriers[["participant_id", "age", "age_of_onset"]]
info = pd.merge(info, age_info, 
                       on="participant_id", how="left")

# Participants with more than one row per gene per ICD-10 code (should not occur)
duplicates = info[info.duplicated(['participant_id', 'gene', 'icd10'], keep = False)]

# Dropping these duplicates
info = info.drop_duplicates(['participant_id', 'gene', 'icd10'])

# Dropping any duplicate rows
info = info.drop_duplicates()

# Cleaning age columns
info.loc[info['age_of_onset']==-999, 'age_of_onset'] = 0
info.loc[info['age_of_onset']==999, 'age_of_onset'] = np.nan

# Adding column that specifies whether the participant is below the typical age of onset for the disease (1=yes)
# create a list of conditions (age of onset for each disease)
conditions = [
    (info['variant_specific_disease'] == 'Hypertrophic Cardiomyopathy') & (info['age'] <= 25),
    (info['variant_specific_disease'] == 'Dilated Cardiomyopathy'),
    (info['variant_specific_disease'] == 'Dilated Cardiomyopathy / Left Ventricular Noncompaction'),
    (info['variant_specific_disease'] == 'Fabry Disease (Hypertrophic Cardiomyopathy)') & (info['age'] <= 16),
    (info['variant_specific_disease'] == 'Arrhythmogenic Right Ventricular Cardiomyopathy') & (info['age'] <= 26),
    (info['variant_specific_disease'] == 'Long QT / Brugada Syndrome')
]

# create a list of the diseases to assign for each condition
values = [1, 0, 0, 1, 1, 0] # DCM and LQTS all ages are onset
# create a new column and use np.select to assign values to it using our lists as arguments
info['below_onset'] = np.select(conditions, values)

# setting below onset column to 0 if age of onset is specified
info.loc[info['age_of_onset'].notnull(), 'below_onset'] = 0

# Some discrepancy in some rows, fixing manually
info.loc[(info['variant_specific_disease'] == 'Long QT / Brugada Syndrome'), 'below_onset'] = 0
info.loc[(info['variant_specific_disease'] == 'Dilated Cardiomyopathy'), 'below_onset'] = 0
info.loc[(info['variant_specific_disease'] == 'Dilated Cardiomyopathy / Left Ventricular Noncompaction'), 'below_onset'] = 0

# Dropping duplicates
info = info.drop_duplicates()

# Checking to see if there are any rare disease relatives that are recruited for the variant disease (proband was recruited for that disease), and have no related ICD10 codes
# If this was true for a participant, cannot guarantee they display the disease if penetrance is based on recruitment disease
info.loc[(info['participant_id']=='Rare Disease Relative') &
                    (info['recruited_for_variant_disease']==1.0) & 
                    (info['variant_related_icd10']==0)] # no participants meet criteria

# Writing to file
info.to_csv("final_penetrance.txt", index = False)

# Splitting out into diseases
age_info = info[['participant_id', 'age', 'gene']]
age_info = age_info.drop_duplicates()

HCM = age_info.loc[(age_info['gene'] == 'MYBPC3') | 
    (age_info['gene'] == 'MYH7') |
    (age_info['gene'] == 'TNNI3') |
    (age_info['gene'] == 'TPM1') |
    (age_info['gene'] == 'MYL3') |
    (age_info['gene'] == 'ACTC1') |
    (age_info['gene'] == 'MYL2') |
    (age_info['gene'] == 'PRKAG2')]
DCM = age_info.loc[(age_info['gene'] == 'TNNT2') | 
                            (age_info['gene'] == 'LMNA')]
FD = age_info.loc[age_info['gene'] == 'GLA']
ARVC = age_info.loc[(age_info['gene'] == 'PKP2') |
                             (age_info['gene'] == 'DSP') |
                             (age_info['gene'] == 'DSC2') |
                             (age_info['gene'] == 'TMEM43') |
                             (age_info['gene'] == 'DSG2')]
LQTS_BS = age_info.loc[(age_info['gene'] == 'KCNQ1') |
                                (age_info['gene'] == 'KCNH2') |
                                (age_info['gene'] == 'SCN5A')]

HCM = HCM.sort_values('participant_id')
DCM = DCM.sort_values('participant_id')
FD = FD.sort_values('participant_id')
ARVC = ARVC.sort_values('participant_id')
LQTS_BS = LQTS_BS.sort_values('participant_id')

# Creating new disease dictionary with dataframes
phenotype_dict = {'HCM': HCM,
                 'DCM': DCM,
                 'FD': FD,
                 'ARVC': ARVC,
                 'LQTS_BS':LQTS_BS}

# Looking at total carriers, average age, and age IQR
for key in phenotype_dict:
    print(key, ':')
    print('')
    
    average_age = phenotype_dict[key]['age'].mean()
    q3 = np.percentile(phenotype_dict[key]['age'], [75])
    q1 = np.percentile(phenotype_dict[key]['age'], [25])
    
    print('Total carriers:', len(phenotype_dict[key]['participant_id'].unique()))
    print('Average age: ', average_age)
    print('IQR: ', q1, '-', q3)
    print('')
    print('-----------------------------------')
    print('')

# Age info by gene
# Subsetting dataframes for each gene of interest and cleaning

MYBPC3 = age_info.loc[age_info['gene'] == 'MYBPC3']
MYH7 = age_info.loc[age_info['gene'] == 'MYH7']
TNNT2 = age_info.loc[age_info['gene'] == 'TNNT2']
TNNI3 = age_info.loc[age_info['gene'] == 'TNNI3']
TPM1 = age_info.loc[age_info['gene'] == 'TPM1'] 
MYL3 = age_info.loc[age_info['gene'] == 'MYL3'] 
ACTC1 = age_info.loc[age_info['gene'] == 'ACTC1']
PRKAG2 = age_info.loc[age_info['gene'] == 'PRKAG2']
GLA = age_info.loc[age_info['gene'] == 'GLA']
MYL2 = age_info.loc[age_info['gene'] == 'MYL2']
LMNA = age_info.loc[age_info['gene'] == 'LMNA']
PKP2 = age_info.loc[age_info['gene'] == 'PKP2']
DSP = age_info.loc[age_info['gene'] == 'DSP']
DSC2 = age_info.loc[age_info['gene'] == 'DSC2']
TMEM43 = age_info.loc[age_info['gene'] == 'TMEM43']
DSG2 = age_info.loc[age_info['gene'] == 'DSG2']
KCNQ1 = age_info.loc[age_info['gene'] == 'KCNQ1']
KCNH2 = age_info.loc[age_info['gene'] == 'KCNH2']
SCN5A = age_info.loc[age_info['gene'] == 'SCN5A']

MYBPC3 = MYBPC3.sort_values('participant_id')
MYH7 = MYH7.sort_values('participant_id')
TNNT2 = TNNT2.sort_values('participant_id')
TNNI3 = TNNI3.sort_values('participant_id')
TPM1 = TPM1.sort_values('participant_id')
MYL3 = MYL3.sort_values('participant_id')
ACTC1 = ACTC1.sort_values('participant_id')
PRKAG2 = PRKAG2.sort_values('participant_id')
GLA = GLA.sort_values('participant_id')
MYL2 = MYL2.sort_values('participant_id')
LMNA = LMNA.sort_values('participant_id')
PKP2 = PKP2.sort_values('participant_id')
DSP = DSP.sort_values('participant_id')
DSC2 = DSC2.sort_values('participant_id')
TMEM43 = TMEM43.sort_values('participant_id')
DSG2 = DSG2.sort_values('participant_id')
KCNQ1 = KCNQ1.sort_values('participant_id')
KCNH2 = KCNH2.sort_values('participant_id')
SCN5A = SCN5A.sort_values('participant_id')

# Creating dictionary with each gene dataframe
gene_dict = {'MYBPC3': MYBPC3,
         'MYH7': MYH7,
         'TNNT2': TNNT2,
         'TNNI3': TNNI3,
         'TPM1': TPM1, 
         'MYL3': MYL3, 
         'ACTC1': ACTC1,
         'PRKAG2': PRKAG2,
         'GLA': GLA,
         'MYL2': MYL2,
         'LMNA': LMNA,
         'PKP2': PKP2,
         'DSP': DSP,
         'DSC2': DSC2,
         'TMEM43': TMEM43,
         'DSG2': DSG2,
         'KCNQ1': KCNQ1,
         'KCNH2': KCNH2,
         'SCN5A': SCN5A}

# Looking at total carriers, average age and age IQR for each gene
for key in gene_dict:
    print(key, ':')
    print('')
    
    total = len(gene_dict[key]['participant_id'].unique())
    print('Total carriers:', total)
    
    if total > 0:
        average_age = gene_dict[key]['age'].mean()
        q3 = np.percentile(gene_dict[key]['age'], [75])
        q1 = np.percentile(gene_dict[key]['age'], [25])
        
        print('Average age: ', average_age)
        print('IQR: ', q1, '-', q3)
    
    else:
        pass
    
    print('')
    print('-----------------------------------')
    print('')

# creating penetrance table
data = info.groupby(['participant_id', 
                                 'participant_type', 
                                 'gene']).recruited_for_variant_disease.sum().reset_index()

data2 = info.groupby(['participant_id', 
                                 'participant_type', 
                                 'gene']).variant_related_icd10.sum().reset_index()
data2 = data2.drop(columns = ['participant_type', 'gene'])

data3 = info.groupby(['participant_id', 
                                 'participant_type', 
                                 'gene']).below_onset.sum().reset_index()
data3 = data3.drop(columns = ['participant_type', 'gene'])


pen = pd.merge(data, data2, how='left', on='participant_id')
penetrance_indicators = pd.merge(pen, data3, how='left', on='participant_id')
penetrance_indicators = penetrance_indicators.drop_duplicates()
penetrance_indicators.head()

# Participants with more than one row 
duplicates = penetrance_indicators[penetrance_indicators.duplicated(['participant_id'], keep = False)]
duplicates.head(30)

# Grouping by ID, type, and gene to ensure 1 row per participant
# summing indicator columns per participant (any sum above 1 will mean yes)
df = penetrance_indicators.groupby(['participant_id', 
                                    'participant_type', 
                                    'gene'])['recruited_for_variant_disease', 
                                             'variant_related_icd10', 
                                             'below_onset'].sum().reset_index()

# Participants with more than one row now (more than one variant carried)
duplicates = df[df.duplicated(['participant_id'], keep = False)]
duplicates.head(30)

# setting indicator columns to 1 if above 1
df.loc[df['recruited_for_variant_disease']>1.0, 'recruited_for_variant_disease'] = 1.0
df.loc[df['variant_related_icd10']>1.0, 'variant_related_icd10'] = 1.0
df.loc[df['below_onset']>1.0, 'below_onset'] = 1.0
df.head(10)

# Adding columns to phenotype DFs that indicates if participant should be excluded from penetrance analysis
# (below age of onset and have not been recruited and have no related ICD-10s)
# 1 = yes
df['exclude'] = np.nan

df['exclude'] = np.where(((df['recruited_for_variant_disease'] == 0) & 
                         (df['variant_related_icd10'] == 0) & 
                         (df['below_onset'] == 1)), 
                         1, 
                         df['exclude']) 

df['exclude'] = np.where((df['exclude'].isnull()), 
                         0, 
                         df['exclude']) 

# Adding column that indicates if participant was not recruited for the variant disease but has related
# ICD-10 codes (1=yes)
df['only_icd10'] = np.nan
df['only_icd10'] = np.where((df['recruited_for_variant_disease'] == 0) & 
                                    (df['variant_related_icd10'] == 1), 
                                    1, 
                                    df['only_icd10']) 
df['only_icd10'] = np.where((df['only_icd10'].isnull()), 
                                    0, 
                                    df['only_icd10']) 

# Genes with no age of onset - setting below onset to 0 (no)
df.loc[df['gene']=='TNNT2', 'below_onset'] = 0.0
df.loc[df['gene']=='LMNA', 'below_onset'] = 0.0
df.loc[df['gene']=='KCNQ1', 'below_onset'] = 0.0
df.loc[df['gene']=='KCNH2', 'below_onset'] = 0.0
df.loc[df['gene']=='SCN5A', 'below_onset'] = 0.0

df.loc[df['gene']=='TNNT2', 'exclude'] = 0.0
df.loc[df['gene']=='LMNA', 'exclude'] = 0.0
df.loc[df['gene']=='KCNQ1', 'exclude'] = 0.0
df.loc[df['gene']=='KCNH2', 'exclude'] = 0.0
df.loc[df['gene']=='SCN5A', 'exclude'] = 0.0

# Dropping duplicate rows
indicators = df.drop_duplicates()

# Writing to file
indicators.to_csv("carrier_indicators.txt", index = False)

# Splitting out into diseases
HCM = indicators.loc[(indicators['gene'] == 'MYBPC3') | 
    (indicators['gene'] == 'MYH7') |
    (indicators['gene'] == 'TNNI3') |
    (indicators['gene'] == 'TPM1') |
    (indicators['gene'] == 'MYL3') |
    (indicators['gene'] == 'ACTC1') |
    (indicators['gene'] == 'MYL2') |
    (indicators['gene'] == 'PRKAG2')]
DCM = indicators.loc[(indicators['gene'] == 'TNNT2') | 
                            (indicators['gene'] == 'LMNA')]
FD = indicators.loc[indicators['gene'] == 'GLA']
ARVC = indicators.loc[(indicators['gene'] == 'PKP2') |
                             (indicators['gene'] == 'DSP') |
                             (indicators['gene'] == 'DSC2') |
                             (indicators['gene'] == 'TMEM43') |
                             (indicators['gene'] == 'DSG2')]
LQTS_BS = indicators.loc[(indicators['gene'] == 'KCNQ1') |
                                (indicators['gene'] == 'KCNH2') |
                                (indicators['gene'] == 'SCN5A')]

HCM = HCM.sort_values('participant_id')
DCM = DCM.sort_values('participant_id')
FD = FD.sort_values('participant_id')
ARVC = ARVC.sort_values('participant_id')
LQTS_BS = LQTS_BS.sort_values('participant_id')

# Creating new dictionary with dataframes
phenotype_dict = {'HCM': HCM,
                 'DCM': DCM,
                 'FD': FD,
                 'ARVC': ARVC,
                 'LQTS_BS':LQTS_BS}

# ## Assessing penetrance
# Genes of interest
genes = ['MYBPC3', 'MYH7', 'TNNT2', 'TNNI3', 'TPM1', 'MYL3', 'ACTC1', 'PRKAG2', 'GLA', 'MYL2', 'LMNA', 
         'PKP2', 'DSP', 'DSC2', 'TMEM43', 'DSG2', 'KCNQ1', 'KCNH2', 'SCN5A']

# For each disease:
# total carriers
# total carriers recruited for the variant disease
# total carriers with related ICD-10 codes, only ICD-10 codes
# total carriers recruited for AND/OR have related ICD-10 codes

for key in phenotype_dict:
    
    print('Total Number of ', key, 'carriers: ', len(phenotype_dict[key]['participant_id'].unique()))
    
    print('')
    
    print('Number of carriers recruited for ', key, '-related disease:')
    print(len(phenotype_dict[key].loc[phenotype_dict[key]['recruited_for_variant_disease'] == 1.0]['participant_id'].unique()))
    
    print('')
    
    print('Number of carriers with ', key, '-related ICD-10 codes:')
    print(len(phenotype_dict[key].loc[phenotype_dict[key]['variant_related_icd10'] == 1]['participant_id'].unique()))
    
    print('')
    
    print('Number of carriers ONLY recruited for ', key, '-related disease:')
    print(len(phenotype_dict[key].loc[(phenotype_dict[key]['recruited_for_variant_disease'] == 1.0) & 
                                      (phenotype_dict[key]['variant_related_icd10'] == 0)]['participant_id'].unique()))
    
    print('')
    
    print('Number of carriers with ONLY ', key, '-related ICD-10 codes:')
    print(len(phenotype_dict[key].loc[phenotype_dict[key]['only_icd10'] == 1.0]))
    
    print('')
    
    print('Number either recruited for AND/OR have ', key, '-related  ICD-10 codes:')
    print(len(phenotype_dict[key].loc[((phenotype_dict[key]['recruited_for_variant_disease'] == 1.0) & (phenotype_dict[key]['variant_related_icd10'] == 0)) | 
                                      ((phenotype_dict[key]['recruited_for_variant_disease'] == 0.0) & (phenotype_dict[key]['variant_related_icd10'] == 1)) |
                                      ((phenotype_dict[key]['recruited_for_variant_disease'] == 1.0) & (phenotype_dict[key]['variant_related_icd10'] == 1))]['participant_id'].unique()))
    
    print('')

    print('')
    print("--------------------------------------------------------------------------------------------")

# For each gene:
for i in genes:
    
    print('Total Number of ', i, 'carriers: ', len(indicators.loc[(indicators['gene']==i)]['participant_id'].unique()))
    
    print("")
    
    print('Number of carriers recruited for ', i, '-related disease:')
    print(len(indicators.loc[(indicators['gene'] == i) & 
                       (indicators['recruited_for_variant_disease'] == 1.0)]['participant_id'].unique()))
    
    print("")
    
    print('Number of carriers with ', i, '-related ICD-10 codes:')
    print(len(indicators.loc[(indicators['gene'] == i) & 
                       (indicators['variant_related_icd10'] == 1)]['participant_id'].unique()))
    
    print("")
    
    print('Number of carriers ONLY recruited for ', i, '-related disease:')
    print(len(indicators.loc[(indicators['gene'] == i) & 
                       (indicators['recruited_for_variant_disease'] == 1.0) & 
                             (indicators['variant_related_icd10'] == 0)]['participant_id'].unique()))
    
    print("")
    
    print('Number of carriers with ONLY ', i, '-related ICD-10 codes:')
    print(len(indicators.loc[(indicators['gene'] == i) & (indicators['only_icd10'] == 1.0)]))
    
    print("")
    
    print('Number either recruited for AND/OR have ', i, '-related  ICD-10 codes:')
    print(len(indicators.loc[(indicators['gene'] == i) & 
                      (((indicators['recruited_for_variant_disease'] == 1.0) & (indicators['variant_related_icd10'] == 0)) |
                       ((indicators['recruited_for_variant_disease'] == 0.0) & (indicators['variant_related_icd10'] == 1)) |
                       ((indicators['recruited_for_variant_disease'] == 1.0) & (indicators['variant_related_icd10'] == 1)))]['participant_id'].unique()))
    

    print("")
    print("--------------------------------------------------------------------------------------------")

# For each disease:
# total carriers
# total carriers below age of onset
# total carriers below age of onset and not recruited and no ICD-10 codes (to exclude)
# penetrance based on recruitment
# penetrance based on recruitment and/or ICD-10 codes

for key in phenotype_dict:
    
    total = len(phenotype_dict[key]['participant_id'].unique())
    below = len(phenotype_dict[key].loc[phenotype_dict[key]['below_onset'] == 1]['participant_id'].unique())
    below_exclude = len(phenotype_dict[key].loc[phenotype_dict[key]['exclude'] == 1]['participant_id'].unique())
    
    print('')
    print(key, ':')
    print('')
    print('Total number of carriers: ', total)
    print('Number of participants below average age of onset: ', below)
    print('Number of participants below average age of onset not recruited & no related ICD-10s: ', below_exclude)
    print('')
    
    if total == 0:
        
        print('')
        print("--------------------------------------------------------------------------------------------")
        pass
    
    elif total >= 1 and below == 0:
        
        print('')
        
        recruited_carriers = len(phenotype_dict[key].loc[phenotype_dict[key]['recruited_for_variant_disease'] == 1.0]['participant_id'].unique())
        
        print(key, '- Penetrance (based on recruitment):', (recruited_carriers/total))
        
        print('')
        
        and_or_carriers = len(phenotype_dict[key].loc[((phenotype_dict[key]['recruited_for_variant_disease'] == 1.0) & (phenotype_dict[key]['variant_related_icd10'] == 0)) |
                                              ((phenotype_dict[key]['recruited_for_variant_disease'] == 0.0) & (phenotype_dict[key]['variant_related_icd10'] == 1)) |
                                              ((phenotype_dict[key]['recruited_for_variant_disease'] == 1.0) & (phenotype_dict[key]['variant_related_icd10'] == 1))]['participant_id'].unique())
        print(i, '- Penetrance (based on recruitment and/or ICD-10 codes):', (and_or_carriers/total))
        
        print('')
        print("--------------------------------------------------------------------------------------------")
    
    elif total>= 1 and below >= 1:
        
        print('')
        recruited_no_exclude = len(phenotype_dict[key].loc[(phenotype_dict[key]['exclude'] == 0) &
                                                  (phenotype_dict[key]['recruited_for_variant_disease'] == 1.0)]['participant_id'].unique())
        print(key, '- Penetrance (based on recruitment, excluding those below average age of onset):', (recruited_no_exclude/(total-below_exclude)))
        print('')
        
        and_or_exclude = len(phenotype_dict[key].loc[(phenotype_dict[key]['exclude'] == 0) &
                                            (((phenotype_dict[key]['recruited_for_variant_disease'] == 1.0) & (phenotype_dict[key]['variant_related_icd10'] == 0)) |
                                             ((phenotype_dict[key]['recruited_for_variant_disease'] == 0.0) & (phenotype_dict[key]['variant_related_icd10'] == 1)) |
                                             ((phenotype_dict[key]['recruited_for_variant_disease'] == 1.0) & (phenotype_dict[key]['variant_related_icd10'] == 1)))]['participant_id'].unique())
        print(key, '- Penetrance (based on recruitment and/or ICD-10, excluding those below average age of onset):', (and_or_exclude/(total-below_exclude)))
        
        print('')
        print("--------------------------------------------------------------------------------------------")
        print('')
            
    else:
        pass

# For each gene:

for i in genes:
    
    total = len(indicators.loc[(indicators['gene']==i)]['participant_id'].unique())
    below = len(indicators.loc[(indicators['gene'] == i) & 
                               (indicators['below_onset'] == 1)]['participant_id'].unique())
    below_exclude = len(indicators.loc[(indicators['gene'] == i) & 
                                       (indicators['exclude'] == 1)]['participant_id'].unique())
    print('')
    print(i, ':')
    print('')
    print('Total number of carriers: ', total)
    print('Number of participants below average age of onset: ', below)
    print('Number of participants below average age of onset not recruited & no related ICD-10s: ', below_exclude)
    print('')
    
    if total == 0:
        
        print('')
        print("--------------------------------------------------------------------------------------------")
        pass
    
    elif total >= 1 and below == 0:
        
        print('')
        
        recruited_carriers = len(indicators.loc[(indicators['gene'] == i) & 
                                                (indicators['recruited_for_variant_disease'] == 1.0)]['participant_id'].unique())
        
        print(i, '- Penetrance (based on recruitment):', (recruited_carriers/total))
        
        print('')
        
        and_or_carriers = len(indicators.loc[(indicators['gene'] == i) & 
                                             (((indicators['recruited_for_variant_disease'] == 1.0) & (indicators['variant_related_icd10'] == 0)) |
                                              ((indicators['recruited_for_variant_disease'] == 0.0) & (indicators['variant_related_icd10'] == 1)) |
                                              ((indicators['recruited_for_variant_disease'] == 1.0) & (indicators['variant_related_icd10'] == 1)))]['participant_id'].unique())
        print(i, '- Penetrance (based on recruitment and/or ICD-10 codes):', (and_or_carriers/total))
        
        print('')
        print("--------------------------------------------------------------------------------------------")
    
    elif total>= 1 and below >= 1:
        
        print('')
        recruited_no_exclude = len(indicators.loc[(indicators['gene'] == i) & 
                                                  (indicators['exclude'] == 0) &
                                                  (indicators['recruited_for_variant_disease'] == 1.0)]['participant_id'].unique())
        print(i, '- Penetrance (based on recruitment, excluding those below average age of onset):', (recruited_no_exclude/(total-below_exclude)))
        print('')
        
        and_or_exclude = len(indicators.loc[(indicators['gene'] == i) & 
                                            (indicators['exclude'] == 0) &
                                            (((indicators['recruited_for_variant_disease'] == 1.0) & (indicators['variant_related_icd10'] == 0)) |
                                             ((indicators['recruited_for_variant_disease'] == 0.0) & (indicators['variant_related_icd10'] == 1)) |
                                             ((indicators['recruited_for_variant_disease'] == 1.0) & (indicators['variant_related_icd10'] == 1)))]['participant_id'].unique())
        print(i, '- Penetrance (based on recruitment and/or ICD-10, excluding those below average age of onset):', (and_or_exclude/(total-below_exclude)))
        
        print('')
        print("--------------------------------------------------------------------------------------------")
        print('')
            
    else:
        
        pass

# # Visualisation

# Number of carriers per gene
counts = indicators['gene'].value_counts()

# Subsetting wanted information
data = df[['gene', 
           'recruited_for_variant_disease', 
           'variant_related_icd10',
           'only_icd10',
           'below_onset']]

# Grouping by gene to identify counts per indicator
data_2 = data.groupby(['gene'])['recruited_for_variant_disease', 
                                'variant_related_icd10',
                                'only_icd10',
                                'below_onset'].sum().reset_index()

# Creating dataframe for penetrance counts, finding total carriers per gene first
totals = counts.to_frame()
tot = totals.reset_index()
tot = tot.rename(columns={"index": "gene", "gene": "total_carriers"})

# Merging total counts with penetrance counts
counts_data = pd.merge(data_2, tot, how='left', on='gene')
counts_data = counts_data.set_index('gene')

# Creating and/or column and calculating penetrance estimates
counts_data['rec_and_or_icd10'] = counts_data['recruited_for_variant_disease'] + counts_data['only_icd10']
counts_data['recruitment_penetrance'] = counts_data['recruited_for_variant_disease'] / (counts_data['total_carriers'] - counts_data['below_onset'])
counts_data['rec_icd10_penetrance'] = counts_data['rec_and_or_icd10'] / (counts_data['total_carriers'] - counts_data['below_onset'])
counts_data = counts_data.reset_index()
g_penetrance = counts_data[['gene', 'recruitment_penetrance', 'rec_icd10_penetrance']]
g_count = counts_data[['recruited_for_variant_disease', 'variant_related_icd10', 'only_icd10', 'below_onset', 'total_carriers']]

# Cleaning table, subsetting wanted info
g_counts_cleaned = g_count.apply(lambda x: [y if y >= 5 else '<5' for y in x])
gene_pen_table = pd.concat([g_penetrance, g_counts_cleaned], axis=1)
gene_pen_table = gene_pen_table[['gene', 'recruited_for_variant_disease', 'variant_related_icd10', 
                                       'only_icd10', 'below_onset', 'total_carriers', 'recruitment_penetrance', 'rec_icd10_penetrance']]

# Writing to file
gene_pen_table.to_csv("gene_penetrance_table.txt", index = False)

# Creating column indicating what disease each gene belongs to
disease_counts_data = counts_data
disease_counts_data['variant_disease'] = np.where(((disease_counts_data['gene'] == 'MYBPC3') | 
                                    (disease_counts_data['gene'] == 'MYH7') | 
                                    (disease_counts_data['gene'] == 'TNNI3') | 
                                    (disease_counts_data['gene'] == 'TPM1') | 
                                    (disease_counts_data['gene'] == 'MYL3') | 
                                    (disease_counts_data['gene'] == 'ACTC1') | 
                                    (disease_counts_data['gene'] == 'MYL2') | 
                                    (disease_counts_data['gene'] == 'PRKAG2')) , 'HCM', np.nan)
disease_counts_data['variant_disease'] = np.where(((disease_counts_data['gene'] == 'TNNT2') | 
                                    (disease_counts_data['gene'] == 'LMNA')) , 'DCM', disease_counts_data['variant_disease'])
disease_counts_data['variant_disease'] = np.where((disease_counts_data['gene'] == 'GLA'), 'FD', disease_counts_data['variant_disease'])
disease_counts_data['variant_disease'] = np.where(((disease_counts_data['gene'] == 'PKP2') |
                                    (disease_counts_data['gene'] == 'DSP') |
                                    (disease_counts_data['gene'] == 'DSC2') |
                                    (disease_counts_data['gene'] == 'TMEM43') |
                                    (disease_counts_data['gene'] == 'DSG2')), 'ARVC', disease_counts_data['variant_disease'])
disease_counts_data['variant_disease'] = np.where(((disease_counts_data['gene'] == 'KCNQ1') |
                                    (disease_counts_data['gene'] == 'KCNH2') |
                                    (disease_counts_data['gene'] == 'SCN5A')), 'LQTS', disease_counts_data['variant_disease'])

# Grouping by disease
disease_counts_data = disease_counts_data.drop(columns= ['gene', 'recruitment_penetrance', 'rec_icd10_penetrance'])
disease_counts_data = disease_counts_data.groupby('variant_disease').sum()
disease_counts_data = disease_counts_data.reset_index()

# Creating and/or column and calculating penetrance estimates for diseases
disease_counts_data['rec_and_or_icd10'] = disease_counts_data['recruited_for_variant_disease'] + disease_counts_data['only_icd10']
disease_counts_data['recruitment_penetrance'] = disease_counts_data['recruited_for_variant_disease'] / (disease_counts_data['total_carriers'] - disease_counts_data['below_onset'])
disease_counts_data['rec_icd10_penetrance'] = disease_counts_data['rec_and_or_icd10'] / (disease_counts_data['total_carriers'] - disease_counts_data['below_onset'])
penetrance = disease_counts_data[['variant_disease', 'recruitment_penetrance', 'rec_icd10_penetrance']]
disease_count = disease_counts_data[['recruited_for_variant_disease', 'variant_related_icd10', 'only_icd10', 'below_onset', 'total_carriers']]

# Cleaning table
disease_counts_cleaned = disease_count.apply(lambda x: [y if y >= 5 else '<5' for y in x])

# Adding total carrier info to penetrance info, subsetting wanted info
disease_pen_table = pd.concat([penetrance, disease_counts_cleaned], axis=1)
disease_pen_table = disease_pen_table[['variant_disease', 'recruited_for_variant_disease', 'variant_related_icd10', 
                                       'only_icd10', 'below_onset', 'total_carriers', 'recruitment_penetrance', 'rec_icd10_penetrance']]

# Writing to file
disease_pen_table.to_csv("disease_penetrance_table.txt", index = False)