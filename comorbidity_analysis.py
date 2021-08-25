import labkey
import pandas as pd
import numpy as np

# COMORBIDITY ANALYSIS

# File with all identified comorbidities for each participant
diseases = pd.read_csv('308_affected.txt.gz', sep='\t')

# Reading in file with carrier info
carrier_status = pd.read_csv('full_carriers_info.txt', sep=',')

# Genes of interest
genes = ['MYBPC3', 'MYH7', 'TNNT2', 'TNNI3', 'TPM1', 'MYL3', 'ACTC1', 'PRKAG2', 'GLA', 'MYL2', 'LMNA', 
         'PKP2', 'DSP', 'DSC2', 'TMEM43', 'DSG2', 'KCNQ1', 'KCNH2', 'SCN5A']

# Subsetting wanted information
carriers = carrier_status[['participant_id', 'gene', 'participant_type', 'age']]
carriers = carriers.loc[carriers['gene'].isin(genes)]

# Reading in all participants with genomic info and subsetting to wanted info
participants = pd.read_csv("participants_in_gvcf.txt", sep = "\t")
all_info = participants[['participant_id', 'participant_type', 'participant_phenotypic_sex', 
                         'participant_ethnic_category_simple', 'year_of_birth', 'date_of_death', 'related']]

# Adding death info
all_info['died'] = np.where(pd.notnull(all_info['date_of_death']), 1, None)

# Extracting year of death
all_info['year_of_death'] = pd.DatetimeIndex(all_info['date_of_death']).year

# Calculate age
all_info['age'] = ((np.where(pd.isnull(all_info['date_of_death']), 2021, all_info['year_of_death']) - all_info['year_of_birth']))
all_info = all_info

# Subsetting wanted info
participant_info =  all_info[['participant_id', 'participant_type', 'participant_phenotypic_sex', 
                              'participant_ethnic_category_simple', 'died', 'age', 'related']]

# Merging demographic info with carrier info
data = pd.merge(carriers, participant_info, how='outer', on=['participant_id', 'participant_type', 'age'])
data = data.sort_values('participant_id')
data = data.reset_index()
data = data.drop(columns='index')

# Adding carrier status for all genes
data['carrier_status'] = np.where(data['gene'].isin(genes), 1, 0)

# Adding column that indicates what disease the gene belongs to
data['variant_disease'] = np.where(((data['gene'] == 'MYBPC3') | 
                                    (data['gene'] == 'MYH7') | 
                                    (data['gene'] == 'TNNI3') | 
                                    (data['gene'] == 'TPM1') | 
                                    (data['gene'] == 'MYL3') | 
                                    (data['gene'] == 'ACTC1') | 
                                    (data['gene'] == 'MYL2') | 
                                    (data['gene'] == 'PRKAG2')) , 'HCM', np.nan)
data['variant_disease'] = np.where(((data['gene'] == 'TNNT2') | 
                                    (data['gene'] == 'LMNA')) , 'DCM', data['variant_disease'])
data['variant_disease'] = np.where((data['gene'] == 'GLA'), 'FD', data['variant_disease'])
data['variant_disease'] = np.where(((data['gene'] == 'PKP2') |
                                    (data['gene'] == 'DSP') |
                                    (data['gene'] == 'DSC2') |
                                    (data['gene'] == 'TMEM43') |
                                    (data['gene'] == 'DSG2')), 'ARVC', data['variant_disease'])
data['variant_disease'] = np.where(((data['gene'] == 'KCNQ1') |
                                    (data['gene'] == 'KCNH2') |
                                    (data['gene'] == 'SCN5A')), 'LQTS', data['variant_disease'])

# Filling null deaths with 0 (no)
data["died"] = data["died"].fillna(0)

# Filtering data so only unrelated participants are used
data_unrelated = data.loc[data['related'] == 'Unrelated']

# Looking at duplicates
duplicates = data_unrelated[data_unrelated.duplicated(['participant_id', 'gene'], keep = False)]

# Removing all duplicate rows where patient has multiple entries
clean_data = data_unrelated.drop_duplicates(['participant_id', 'gene'])
clean_data = clean_data.reset_index()
clean_data = clean_data.drop(columns=['index'])

# Observing if any duplicates
duplicates = clean_data[clean_data.duplicated(keep = False)]

# Merging unrelated demographic information with comorbidity data
comorbidities = pd.merge(clean_data, diseases, on='participant_id', how='left')
comorbidities = comorbidities[['participant_id', 'gene', 'participant_type', 'participant_phenotypic_sex', 
                               'participant_ethnic_category_simple', 'age', 
                               'disease', 'category', 'date', 'n_obs', 'carrier_status' , 'related', 'variant_disease']]

# Looking at duplicates
duplicates = comorbidities[comorbidities.duplicated(['participant_id', 'gene', 'disease'], keep = False)]

# Looking at top comorbidities per gene with more than 10 carriers
for i in genes:
    carriers = comorbidities.loc[(comorbidities['carrier_status'] == 1) & (comorbidities['gene'] == i)]
    carrier_total = len(carriers.participant_id.unique())
    
    comorbid_counts = pd.DataFrame({'count': carriers.groupby(['disease']).size()}).reset_index()
    comorbid_counts = comorbid_counts.sort_values('count',ascending=False).reset_index()
    comorbid_counts = comorbid_counts.drop(columns=['index'])
    
    comorbid_counts['proportion'] = np.nan
    comorbid_counts['proportion'] = (comorbid_counts['count'] / carrier_total).round(2)
    
    print('Top comorbidities (single entry) for', i, ' carriers:')
    print('')
    print(comorbid_counts.loc[comorbid_counts['count']>=10])
         
    print('')
    print('--------------------------------------------------------------------------------------')
    print('')

# All comorbidities with at least 10 carriers in at least 1 gene
['Hypertrophic Cardiomyopathy', 'Hypertension', 'Other or unspecified infectious organisms', 'Lower Respiratory Tract Infections', 
 'Bacterial Diseases (excl TB)', 'Asthma', 'Obesity', 'Gastritis and duodenitis', 'Atrial fibrillation', 'Infections of Other or unspecified organs', 
 'Osteoarthritis (excl spine)', 'Other anaemias', 'Benign neoplasm of colon, rectum, anus and anal canal', 'Cataract', 
 'Diverticular disease of intestine (acute and chronic)', 'Ear and Upper Respiratory Tract Infections', 'Primary Malignancy_Breast', 
 'Infections of the digestive system', 'Iron deficiency anaemia', 'Diabetes']

# # Statistics
### Looking at significant difference in associated comorbidities between carriers and non-carriers
#### https://towardsdatascience.com/building-a-logistic-regression-in-python-step-by-step-becd4d56c9c8

# Subsetting wanted information
stats_df = comorbidities[['participant_id', 'gene', 'participant_type', 'participant_phenotypic_sex', 
                          'participant_ethnic_category_simple', 'age', 'disease', 'carrier_status', 'variant_disease']]

# All possible comorbidities
diseases = stats_df.disease.unique()

# Creating dataframe where each column is each possible comorbidity
disease_df = pd.DataFrame(columns=diseases)

# Adding comorbidity columns to previous dataframe
classification = pd.concat([stats_df[['participant_id', 'gene', 'variant_disease', 'disease']], disease_df], axis=1)

# Editing each comorbidity column to 1 if the participant has that comorbidity in their data
for col in classification.columns[4:]:
    classification[col] = (classification['disease'] == col).astype(int)

# Dropping original comorbidity column
classification = classification.drop(columns=['disease', np.nan])

# Grouping data by participant to ensure 1 row per participant
df = classification.groupby(['participant_id'])[classification.columns[3:]].apply(lambda x : x.astype(int).sum())
df = df.reset_index()

# Participant IDs list
ids = df[['participant_id']]

# If sum was greater than 1, changing back to 1 as columns should have binary inputs (0= does not have condition, 1=yes)
df_cleaned = df.apply(lambda x: [y if y <= 1 else 1 for y in x]).drop(columns='participant_id')

# Adding participant IDs column back to comorbidity data
clean = pd.concat([ids, df_cleaned], axis=1)

# Demographic info to add back onto comorbidity data
df_simple = stats_df.drop(columns='disease')
df_simple = df_simple.drop_duplicates().reset_index()
df_simple = df_simple.drop(columns='index')

# changing order of wanted data
df_simple = df_simple[['participant_id', 'variant_disease', 'carrier_status', 'gene', 'participant_type',
                       'participant_phenotypic_sex', 'participant_ethnic_category_simple', 'age']]

# Merging demographic data with comorbidity data
final_data = pd.merge(df_simple, clean, how='left', on=['participant_id'])

# Dataframe with comorbidity data and genes only
gene_comorbids = final_data.drop(columns=['variant_disease', 'participant_id', 'participant_type', 'participant_phenotypic_sex', 'participant_ethnic_category_simple', 'age'])

# Grouping comorbidity data by gene for carriers only
gene_df = gene_comorbids.groupby(['carrier_status', 'gene'])[classification.columns[3:]].apply(lambda x : x.astype(int).sum())
gene_df = gene_df.reset_index()
gene_df = gene_df.loc[gene_df['carrier_status']==1]
gene_df = gene_df.drop(columns='carrier_status')

# Creating dataframe with comorbidities for diseases only
disease_comorbids = final_data.drop(columns=['gene', 'carrier_status', 'participant_id', 'participant_type', 'participant_phenotypic_sex', 'participant_ethnic_category_simple', 'age'])

# Grouping comorbidity data by disease of interest
disease_df = disease_comorbids.groupby(['variant_disease'])[classification.columns[3:]].apply(lambda x : x.astype(int).sum())
disease_df = disease_df.reset_index()
disease_df = disease_df.loc[disease_df['variant_disease'] != np.nan]

genes = ['MYBPC3', 'MYH7', 'TNNT2', 'TNNI3', 'TPM1', 'MYL3', 'ACTC1', 'PRKAG2', 'GLA', 'MYL2', 'LMNA', 
         'PKP2', 'DSP', 'DSC2', 'TMEM43', 'DSG2', 'KCNQ1', 'KCNH2', 'SCN5A']
# Number of participants per gene with info for comorbidities
for i in genes:
    total = len(final_data.loc[final_data['gene'] == i]['participant_id'].unique())
    print(i, ': ', total)

# Creating dataframe with total carriers per gene
carrier_counts = pd.DataFrame({'MYBPC3': 89,
                   'MYH7': 62, 
                   'TNNT2': 45, 
                   'TNNI3': 15, 
                   'TPM1': 5, 
                   'MYL3': 0, 
                   'ACTC1': 0, 
                   'PRKAG2': 20, 
                   'GLA': 13, 
                   'MYL2': 5, 
                   'LMNA': 57, 
                   'PKP2': 67, 
                   'DSP': 23, 
                   'DSC2': 22, 
                   'TMEM43': 33, 
                   'DSG2': 38, 
                   'KCNQ1': 121, 
                   'KCNH2': 74, 
                   'SCN5A': 58}, 
index = ['total_carriers'])
carrier_counts = carrier_counts.transpose()
carrier_counts = carrier_counts.reset_index()
carrier_counts = carrier_counts.rename(columns={'index': 'gene'})

# Creating column that indicates what disease each gene belongs to
carrier_counts['variant_disease'] = np.where(((carrier_counts['gene'] == 'MYBPC3') | 
                                    (carrier_counts['gene'] == 'MYH7') | 
                                    (carrier_counts['gene'] == 'TNNI3') | 
                                    (carrier_counts['gene'] == 'TPM1') | 
                                    (carrier_counts['gene'] == 'MYL3') | 
                                    (carrier_counts['gene'] == 'ACTC1') | 
                                    (carrier_counts['gene'] == 'MYL2') | 
                                    (carrier_counts['gene'] == 'PRKAG2')) , 'HCM', np.nan)
carrier_counts['variant_disease'] = np.where(((carrier_counts['gene'] == 'TNNT2') | 
                                    (carrier_counts['gene'] == 'LMNA')) , 'DCM', carrier_counts['variant_disease'])
carrier_counts['variant_disease'] = np.where((carrier_counts['gene'] == 'GLA'), 'FD', carrier_counts['variant_disease'])
carrier_counts['variant_disease'] = np.where(((carrier_counts['gene'] == 'PKP2') |
                                    (carrier_counts['gene'] == 'DSP') |
                                    (carrier_counts['gene'] == 'DSC2') |
                                    (carrier_counts['gene'] == 'TMEM43') |
                                    (carrier_counts['gene'] == 'DSG2')), 'ARVC', carrier_counts['variant_disease'])
carrier_counts['variant_disease'] = np.where(((carrier_counts['gene'] == 'KCNQ1') |
                                    (carrier_counts['gene'] == 'KCNH2') |
                                    (carrier_counts['gene'] == 'SCN5A')), 'LQTS', carrier_counts['variant_disease'])

# Grouping gene counts by disease 
disease_counts = pd.DataFrame({'total_carriers': carrier_counts.groupby(['variant_disease'])['total_carriers'].sum()}).reset_index()

# Merging carriers and disease info onto comorbidity data
gene_co_info = pd.merge(carrier_counts, gene_df, how='left', on='gene')

# Merging disease counts to comorbidity data
disease_co_info = pd.merge(disease_counts, disease_df, how='left', on='variant_disease')

disease_names = disease_co_info['variant_disease']

# Looking at top comorbidities per disease with more than 15 carriers
for i in disease_names:
    
    carriers = comorbidities.loc[(comorbidities['carrier_status'] == 1) & (comorbidities['variant_disease'] == i)]
    carrier_total = len(carriers.participant_id.unique())
       
    comorbid_counts = pd.DataFrame({'count': carriers.groupby(['disease']).size()}).reset_index()
    comorbid_counts = comorbid_counts.sort_values('count',ascending=False).reset_index()
    comorbid_counts = comorbid_counts.drop(columns=['index'])
    
    comorbid_counts['proportion'] = np.nan
    comorbid_counts['proportion'] = (comorbid_counts['count'] / carrier_total).round(2)
    
    print('Top comorbidities (single entry) for', i, ' carriers:')
    print('')
    print(comorbid_counts.loc[comorbid_counts['count']>=15])
          
    print('')
    print('--------------------------------------------------------------------------------------')
    print('')

# Subsetting comorbidity data to conditions with at least 15 carriers in one disease
more_than_15 = disease_co_info[['variant_disease', 'total_carriers','Hypertrophic Cardiomyopathy', 
                                    'Other or unspecified infectious organisms', 
                                    'Hypertension',
                                    'Bacterial Diseases (excl TB)', 
                                    'Lower Respiratory Tract Infections', 
                                    'Infections of Other or unspecified organs', 
                                    'Viral diseases (excl chronic hepatitis/HIV)',
                                    'Obesity', 
                                    'Ear and Upper Respiratory Tract Infections', 
                                    'Asthma', 
                                    'Infections of the digestive system', 
                                    'Atrial fibrillation', 
                                    'Gastritis and duodenitis', 
                                    'Urinary Tract Infections', 
                                    'Gastro-oesophageal reflux disease', 
                                    'Coronary heart disease not otherwise specified', 
                                    'Osteoarthritis (excl spine)', 
                                    'Anxiety disorders', 
                                    'Intellectual disability',
                                    'Other anaemias', 
                                    'Cataract', 
                                    'Iron deficiency anaemia', 
                                    'Depression', 
                                    'Infection of skin and subcutaneous tissues', 
                                    'Benign neoplasm of colon, rectum, anus and anal canal', 
                                    'Diabetes', 
                                    'Epilepsy', 
                                    'Primary Malignancy_Breast', 
                                    'Diverticular disease of intestine (acute and chronic)', 
                                    'Hypo or hyperthyroidism', 
                                    'Septicaemia', 
                                    'Abdominal Hernia', 'Gastritis and duodenitis', 
                                'Secondary Malignancy_Lymph Nodes',  
                                'Heart failure', 
                                'Myocardial infarction', 
                                'Acute Kidney Injury', 
                                'Abdominal Hernia', 
                                'Intervertebral disc disorders']]

# Editing dataframe so all values are float integers
ints = more_than_15.iloc[:, 1:].astype(float)
dis = more_than_15[['variant_disease']]
more15 = pd.concat([dis, ints], axis=1)

# Creating dataframe per disease where all numbers represent the proportion of carriers with that condition
arvc = more15.loc[more15['variant_disease'] == 'ARVC'] 
for col in arvc.columns[2:]:
    arvc[col] = arvc[col].div(252)

dcm = more15.loc[more15['variant_disease'] == 'DCM'] 
for col in dcm.columns[2:]:
    dcm[col] = dcm[col].div(155)
    
fd = more15.loc[more15['variant_disease'] == 'FD'] 
for col in fd.columns[2:]:
    fd[col] = fd[col].div(25)
    
hcm = more15.loc[more15['variant_disease'] == 'HCM'] 
for col in hcm.columns[2:]:
    hcm[col] = hcm[col].div(261)
    
lqts = more15.loc[more15['variant_disease'] == 'LQTS'] 
for col in lqts.columns[2:]:
    lqts[col] = lqts[col].div(357)

# Merging all disease comorbidity proportions into one dataframe
m15 = pd.concat([arvc, dcm, fd, hcm, lqts])

# Setting index to disease
disease_proportion_map = m15.set_index('variant_disease').drop(columns='total_carriers')

# Transposing data so diseases are columns
disease_prop_transp = disease_proportion_map.transpose()

# Constructing proportion heatmap for diseases
import matplotlib.pyplot as plt
import seaborn as sns
from pylab import savefig

fig, ax = plt.subplots(figsize=(12,12))
svm = sns.heatmap(disease_prop_transp, annot=True, ax=ax, cbar_kws={'label': 'Proportion of Carriers'})
figure = svm.get_figure()    
figure.savefig('disease_proportion_heatmap.png', dpi=400)

# Comorbidities with more than 10 carriers in at least 1 gene, subsetting to these conditions
more_than_10 = gene_co_info[['gene', 'total_carriers',
                             'Hypertrophic Cardiomyopathy', 'Hypertension', 'Other or unspecified infectious organisms', 
                             'Lower Respiratory Tract Infections', 'Bacterial Diseases (excl TB)', 'Asthma', 'Obesity', 'Gastritis and duodenitis', 
                             'Atrial fibrillation', 'Infections of Other or unspecified organs', 
                             'Osteoarthritis (excl spine)', 'Other anaemias', 'Benign neoplasm of colon, rectum, anus and anal canal', 'Cataract', 
                             'Diverticular disease of intestine (acute and chronic)', 'Ear and Upper Respiratory Tract Infections', 'Primary Malignancy_Breast', 
                             'Infections of the digestive system', 'Iron deficiency anaemia', 'Diabetes']]

# Calculating proportion of carriers per comorbidity
for col in more_than_10.columns[2:]:
    more_than_10[col] = (more_than_10[col] / more_than_10['total_carriers'])

# Subsetting data to genes with at least 10 carriers
more_than_10_carriers = more_than_10.loc[more_than_10['total_carriers']>10]

# Setting index to gene
proportion_map = more_than_10_carriers.set_index('gene').drop(columns='total_carriers')
proportion_map.head(20)

# Transposing data so genes are columns
prop_transp = proportion_map.transpose()

# Constructing comorbidity proportion heatmap per gene
fig, ax = plt.subplots(figsize=(12,12))
svm = sns.heatmap(prop_transp, annot=True, ax=ax, cbar_kws={'label': 'Proportion of Carriers'})
figure = svm.get_figure()    
figure.savefig('gene_proportion_heatmap.png', dpi=400)

import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

# Setting all decimals to 4 places
pd.options.display.float_format = '{:.4f}'.format
# Setting seed value
from numpy.random import seed
seed(7)

# https://stackoverflow.com/questions/51734180/converting-statsmodels-summary-object-to-pandas-dataframe 

# Function that takes results from statsmodel logistic regression function and constructs a dataframe from it
def results_summary_to_dataframe(results):

    pvals = results.pvalues
    coeff = np.exp(results.params)
    conf_lower = np.exp(results.conf_int()[0])
    conf_higher = np.exp(results.conf_int()[1])

    results_df = pd.DataFrame({"pvals":pvals,
                               "coeff":coeff,
                               "conf_lower":conf_lower,
                               "conf_higher":conf_higher
                                })
    #Reordering...
    results_df = results_df[["coeff","pvals","conf_lower","conf_higher"]]
    return results_df

# Editing data to indicators per gene:
## carrier status indicating if participant is carrying particular gene (0= np, 1= yes)
## participant type: proband = 0, relative = 1, cancer = 2
## gender: 0 = male, 1 = female
## Subsetting into outcome (carrier status) and all covariates (age, gender, participant type, and all comorbidities)
# Building logistic regression and fitting to data
## Subsetting result to those with p-values of 0.05 or lower

# ## MYBPC3
MYBPC3 = final_data
MYBPC3['carrier_status'] = np.where(MYBPC3['gene'] != 'MYBPC3', 0, 1)
MYBPC3 = MYBPC3.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
MYBPC3['participant_type'] = np.where((MYBPC3['participant_type']=='Rare Disease Proband'), 0, 1)
MYBPC3['participant_type'] = np.where((MYBPC3['participant_type']=='Cancer'), 2, MYBPC3['participant_type'])
MYBPC3['participant_phenotypic_sex'] = np.where((MYBPC3['participant_phenotypic_sex']=='Male'), 0, 1)
MYBPC3 = MYBPC3.dropna(how='all')
X = MYBPC3.iloc[:,4:]
y = MYBPC3[['carrier_status']]
 # building the model and fitting the data
log_reg = sm.Logit(y, X).fit(method='bfgs')
results_MYBPC3 = results_summary_to_dataframe(log_reg)
results_MYBPC3 = results_MYBPC3.loc[results_MYBPC3['pvals'] <= 0.05]

# ## MYH7
MYH7 = final_data
MYH7['carrier_status'] = np.where(MYH7['gene'] != 'MYH7', 0, 1)
MYH7 = MYH7.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
MYH7['participant_type'] = np.where((MYH7['participant_type']=='Rare Disease Proband'), 0, 1)
MYH7['participant_type'] = np.where((MYH7['participant_type']=='Cancer'), 2, MYH7['participant_type'])
MYH7['participant_phenotypic_sex'] = np.where((MYH7['participant_phenotypic_sex']=='Male'), 0, 1)
MYH7 = MYH7.dropna(how='all')
X = MYH7.iloc[:,4:]
y = MYH7[['carrier_status']]
 # building the model and fitting the data
log_reg_MYH7 = sm.Logit(y, X).fit(method='bfgs')
results_MYH7 = results_summary_to_dataframe(log_reg_MYH7)
results_MYH7 = results_MYH7.loc[results_MYH7['pvals'] <= 0.05]

# ## TNNT2
TNNT2 = final_data
TNNT2['carrier_status'] = np.where(TNNT2['gene'] != 'TNNT2', 0, 1)
TNNT2 = TNNT2.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
TNNT2['participant_type'] = np.where((TNNT2['participant_type']=='Rare Disease Proband'), 0, 1)
TNNT2['participant_type'] = np.where((TNNT2['participant_type']=='Cancer'), 2, TNNT2['participant_type'])
TNNT2['participant_phenotypic_sex'] = np.where((TNNT2['participant_phenotypic_sex']=='Male'), 0, 1)
TNNT2 = TNNT2.dropna(how='all')
X = TNNT2.iloc[:,4:]
y = TNNT2[['carrier_status']]
 # building the model and fitting the data
log_reg_TNNT2 = sm.Logit(y, X).fit(method='bfgs')
results_TNNT2 = results_summary_to_dataframe(log_reg_TNNT2)
results_TNNT2 = results_TNNT2.loc[results_TNNT2['pvals'] <= 0.05]

# ## TNNI3
TNNI3 = final_data
TNNI3['carrier_status'] = np.where(TNNI3['gene'] != 'TNNI3', 0, 1)
TNNI3 = TNNI3.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
TNNI3['participant_type'] = np.where((TNNI3['participant_type']=='Rare Disease Proband'), 0, 1)
TNNI3['participant_type'] = np.where((TNNI3['participant_type']=='Cancer'), 2, TNNI3['participant_type'])
TNNI3['participant_phenotypic_sex'] = np.where((TNNI3['participant_phenotypic_sex']=='Male'), 0, 1)
TNNI3 = TNNI3.dropna(how='all')
X = TNNI3.iloc[:,4:]
y = TNNI3[['carrier_status']]
 # building the model and fitting the data
log_reg_TNNI3 = sm.Logit(y, X).fit(method='bfgs')
results_TNNI3 = results_summary_to_dataframe(log_reg_TNNI3)
results_TNNI3 = results_TNNI3.loc[results_TNNI3['pvals'] <= 0.05]

# ## TPM1
TPM1 = final_data
TPM1['carrier_status'] = np.where(TPM1['gene'] != 'TPM1', 0, 1)
TPM1 = TPM1.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
TPM1['participant_type'] = np.where((TPM1['participant_type']=='Rare Disease Proband'), 0, 1)
TPM1['participant_type'] = np.where((TPM1['participant_type']=='Cancer'), 2, TPM1['participant_type'])
TPM1['participant_phenotypic_sex'] = np.where((TPM1['participant_phenotypic_sex']=='Male'), 0, 1)
TPM1 = TPM1.dropna(how='all')
X = TPM1.iloc[:,4:]
y = TPM1[['carrier_status']]
 # building the model and fitting the data
log_reg_TPM1 = sm.Logit(y, X).fit(method='bfgs')
results_TPM1 = results_summary_to_dataframe(log_reg_TPM1)
results_TPM1 = results_TPM1.loc[results_TPM1['pvals'] <= 0.05]

# ## MYL3
MYL3 = final_data
MYL3['carrier_status'] = np.where(MYL3['gene'] != 'MYL3', 0, 1)
MYL3 = MYL3.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
MYL3['participant_type'] = np.where((MYL3['participant_type']=='Rare Disease Proband'), 0, 1)
MYL3['participant_type'] = np.where((MYL3['participant_type']=='Cancer'), 2, MYL3['participant_type'])
MYL3['participant_phenotypic_sex'] = np.where((MYL3['participant_phenotypic_sex']=='Male'), 0, 1)
MYL3 = MYL3.dropna(how='all')
X = MYL3.iloc[:,4:]
y = MYL3[['carrier_status']]
 # building the model and fitting the data
log_reg_MYL3 = sm.Logit(y, X).fit(method='bfgs')
results_MYL3 = results_summary_to_dataframe(log_reg_MYL3)
results_MYL3 = results_MYL3.loc[results_MYL3['pvals'] <= 0.05]

# ## ACTC1
ACTC1 = final_data
ACTC1['carrier_status'] = np.where(ACTC1['gene'] != 'ACTC1', 0, 1)
ACTC1 = ACTC1.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
ACTC1['participant_type'] = np.where((ACTC1['participant_type']=='Rare Disease Proband'), 0, 1)
ACTC1['participant_type'] = np.where((ACTC1['participant_type']=='Cancer'), 2, ACTC1['participant_type'])
ACTC1['participant_phenotypic_sex'] = np.where((ACTC1['participant_phenotypic_sex']=='Male'), 0, 1)
ACTC1 = ACTC1.dropna(how='all')
X = ACTC1.iloc[:,4:]
y = ACTC1[['carrier_status']]
 # building the model and fitting the data
log_reg_ACTC1 = sm.Logit(y, X).fit(method='bfgs')
results_ACTC1 = results_summary_to_dataframe(log_reg_ACTC1)
results_ACTC1 = results_ACTC1.loc[results_ACTC1['pvals'] <= 0.05]

# ## PRKAG2
PRKAG2 = final_data
PRKAG2['carrier_status'] = np.where(PRKAG2['gene'] != 'PRKAG2', 0, 1)
PRKAG2 = PRKAG2.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
PRKAG2['participant_type'] = np.where((PRKAG2['participant_type']=='Rare Disease Proband'), 0, 1)
PRKAG2['participant_type'] = np.where((PRKAG2['participant_type']=='Cancer'), 2, PRKAG2['participant_type'])
PRKAG2['participant_phenotypic_sex'] = np.where((PRKAG2['participant_phenotypic_sex']=='Male'), 0, 1)
PRKAG2 = PRKAG2.dropna(how='all')
X = PRKAG2.iloc[:,4:]
y = PRKAG2[['carrier_status']]
 # building the model and fitting the data
log_reg_PRKAG2 = sm.Logit(y, X).fit(method='bfgs')
results_PRKAG2 = results_summary_to_dataframe(log_reg_PRKAG2)
results_PRKAG2 = results_PRKAG2.loc[results_PRKAG2['pvals'] <= 0.05]

# ## GLA
GLA = final_data
GLA['carrier_status'] = np.where(GLA['gene'] != 'GLA', 0, 1)
GLA = GLA.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
GLA['participant_type'] = np.where((GLA['participant_type']=='Rare Disease Proband'), 0, 1)
GLA['participant_type'] = np.where((GLA['participant_type']=='Cancer'), 2, GLA['participant_type'])
GLA['participant_phenotypic_sex'] = np.where((GLA['participant_phenotypic_sex']=='Male'), 0, 1)
GLA = GLA.dropna(how='all')
X = GLA.iloc[:,4:]
y = GLA[['carrier_status']]
 # building the model and fitting the data
log_reg_GLA = sm.Logit(y, X).fit(method='bfgs')
results_GLA = results_summary_to_dataframe(log_reg_GLA)
results_GLA = results_GLA.loc[results_GLA['pvals'] <= 0.05]

# ## MYL2
MYL2 = final_data
MYL2['carrier_status'] = np.where(MYL2['gene'] != 'MYL2', 0, 1)
MYL2 = MYL2.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
MYL2['participant_type'] = np.where((MYL2['participant_type']=='Rare Disease Proband'), 0, 1)
MYL2['participant_type'] = np.where((MYL2['participant_type']=='Cancer'), 2, MYL2['participant_type'])
MYL2['participant_phenotypic_sex'] = np.where((MYL2['participant_phenotypic_sex']=='Male'), 0, 1)
MYL2 = MYL2.dropna(how='all')
X = MYL2.iloc[:,4:]
y = MYL2[['carrier_status']]
 # building the model and fitting the data
log_reg_MYL2 = sm.Logit(y, X).fit(method='bfgs')
results_MYL2 = results_summary_to_dataframe(log_reg_MYL2)
results_MYL2 = results_MYL2.loc[results_MYL2['pvals'] <= 0.05]

# ## LMNA
LMNA = final_data
LMNA['carrier_status'] = np.where(LMNA['gene'] != 'LMNA', 0, 1)
LMNA = LMNA.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
LMNA['participant_type'] = np.where((LMNA['participant_type']=='Rare Disease Proband'), 0, 1)
LMNA['participant_type'] = np.where((LMNA['participant_type']=='Cancer'), 2, LMNA['participant_type'])
LMNA['participant_phenotypic_sex'] = np.where((LMNA['participant_phenotypic_sex']=='Male'), 0, 1)
LMNA = LMNA.dropna(how='all')
X = LMNA.iloc[:,4:]
y = LMNA[['carrier_status']]
 # building the model and fitting the data
log_reg_LMNA = sm.Logit(y, X).fit(method='bfgs')
results_LMNA = results_summary_to_dataframe(log_reg_LMNA)
results_LMNA = results_LMNA.loc[results_LMNA['pvals'] <= 0.05]

# ## PKP2
PKP2 = final_data
PKP2['carrier_status'] = np.where(PKP2['gene'] != 'PKP2', 0, 1)
PKP2 = PKP2.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
PKP2['participant_type'] = np.where((PKP2['participant_type']=='Rare Disease Proband'), 0, 1)
PKP2['participant_type'] = np.where((PKP2['participant_type']=='Cancer'), 2, PKP2['participant_type'])
PKP2['participant_phenotypic_sex'] = np.where((PKP2['participant_phenotypic_sex']=='Male'), 0, 1)
PKP2 = PKP2.dropna(how='all')
X = PKP2.iloc[:,4:]
y = PKP2[['carrier_status']]
 # building the model and fitting the data
log_reg_PKP2 = sm.Logit(y, X).fit(method='bfgs')
results_PKP2 = results_summary_to_dataframe(log_reg_PKP2)
results_PKP2 = results_PKP2.loc[results_PKP2['pvals'] <= 0.05]

# ## DSP
DSP = final_data
DSP['carrier_status'] = np.where(DSP['gene'] != 'DSP', 0, 1)
DSP = DSP.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
DSP['participant_type'] = np.where((DSP['participant_type']=='Rare Disease Proband'), 0, 1)
DSP['participant_type'] = np.where((DSP['participant_type']=='Cancer'), 2, DSP['participant_type'])
DSP['participant_phenotypic_sex'] = np.where((DSP['participant_phenotypic_sex']=='Male'), 0, 1)
DSP = DSP.dropna(how='all')
X = DSP.iloc[:,4:]
y = DSP[['carrier_status']]
 # building the model and fitting the data
log_reg_DSP = sm.Logit(y, X).fit(method='bfgs')
results_DSP = results_summary_to_dataframe(log_reg_DSP)
results_DSP = results_DSP.loc[results_DSP['pvals'] <= 0.05]

# ## DSC2
DSC2 = final_data
DSC2['carrier_status'] = np.where(DSC2['gene'] != 'DSC2', 0, 1)
DSC2 = DSC2.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
DSC2['participant_type'] = np.where((DSC2['participant_type']=='Rare Disease Proband'), 0, 1)
DSC2['participant_type'] = np.where((DSC2['participant_type']=='Cancer'), 2, DSC2['participant_type'])
DSC2['participant_phenotypic_sex'] = np.where((DSC2['participant_phenotypic_sex']=='Male'), 0, 1)
DSC2 = DSC2.dropna(how='all')
X = DSC2.iloc[:,4:]
y = DSC2[['carrier_status']]
 # building the model and fitting the data
log_reg_DSC2 = sm.Logit(y, X).fit(method='bfgs')
results_DSC2 = results_summary_to_dataframe(log_reg_DSC2)
results_DSC2 = results_DSC2.loc[results_DSC2['pvals'] <= 0.05]

# ## TMEM43
TMEM43 = final_data
TMEM43['carrier_status'] = np.where(TMEM43['gene'] != 'TMEM43', 0, 1)
TMEM43 = TMEM43.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
TMEM43['participant_type'] = np.where((TMEM43['participant_type']=='Rare Disease Proband'), 0, 1)
TMEM43['participant_type'] = np.where((TMEM43['participant_type']=='Cancer'), 2, TMEM43['participant_type'])
TMEM43['participant_phenotypic_sex'] = np.where((TMEM43['participant_phenotypic_sex']=='Male'), 0, 1)
TMEM43 = TMEM43.dropna(how='all')
X = TMEM43.iloc[:,4:]
y = TMEM43[['carrier_status']]
 # building the model and fitting the data
log_reg_TMEM43 = sm.Logit(y, X).fit(method='bfgs')
results_TMEM43 = results_summary_to_dataframe(log_reg_TMEM43)
results_TMEM43 = results_TMEM43.loc[results_TMEM43['pvals'] <= 0.05]

# ## DSG2
DSG2 = final_data
DSG2['carrier_status'] = np.where(DSG2['gene'] != 'DSG2', 0, 1)
DSG2 = DSG2.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
DSG2['participant_type'] = np.where((DSG2['participant_type']=='Rare Disease Proband'), 0, 1)
DSG2['participant_type'] = np.where((DSG2['participant_type']=='Cancer'), 2, DSG2['participant_type'])
DSG2['participant_phenotypic_sex'] = np.where((DSG2['participant_phenotypic_sex']=='Male'), 0, 1)
DSG2 = DSG2.dropna(how='all')
X = DSG2.iloc[:,4:]
y = DSG2[['carrier_status']]
 # building the model and fitting the data
log_reg_DSG2 = sm.Logit(y, X).fit(method='bfgs')
results_DSG2 = results_summary_to_dataframe(log_reg_DSG2)
results_DSG2 = results_DSG2.loc[results_DSG2['pvals'] <= 0.05]

# ## KCNQ1
KCNQ1 = final_data
KCNQ1['carrier_status'] = np.where(KCNQ1['gene'] != 'KCNQ1', 0, 1)
KCNQ1 = KCNQ1.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
KCNQ1['participant_type'] = np.where((KCNQ1['participant_type']=='Rare Disease Proband'), 0, 1)
KCNQ1['participant_type'] = np.where((KCNQ1['participant_type']=='Cancer'), 2, KCNQ1['participant_type'])
KCNQ1['participant_phenotypic_sex'] = np.where((KCNQ1['participant_phenotypic_sex']=='Male'), 0, 1)
KCNQ1 = KCNQ1.dropna(how='all')
X = KCNQ1.iloc[:,4:]
y = KCNQ1[['carrier_status']]
 # building the model and fitting the data
log_reg_KCNQ1 = sm.Logit(y, X).fit(method='bfgs')
results_KCNQ1 = results_summary_to_dataframe(log_reg_KCNQ1)
results_KCNQ1 = results_KCNQ1.loc[results_KCNQ1['pvals'] <= 0.05]

# ## KCNH2
KCNH2 = final_data
KCNH2['carrier_status'] = np.where(KCNH2['gene'] != 'KCNH2', 0, 1)
KCNH2 = KCNH2.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
KCNH2['participant_type'] = np.where((KCNH2['participant_type']=='Rare Disease Proband'), 0, 1)
KCNH2['participant_type'] = np.where((KCNH2['participant_type']=='Cancer'), 2, KCNH2['participant_type'])
KCNH2['participant_phenotypic_sex'] = np.where((KCNH2['participant_phenotypic_sex']=='Male'), 0, 1)
KCNH2 = KCNH2.dropna(how='all')
X = KCNH2.iloc[:,4:]
y = KCNH2[['carrier_status']]
 # building the model and fitting the data
log_reg_KCNH2 = sm.Logit(y, X).fit(method='bfgs')
results_KCNH2 = results_summary_to_dataframe(log_reg_KCNH2)
results_KCNH2 = results_KCNH2.loc[results_KCNH2['pvals'] <= 0.05]

# ## SCN5A
SCN5A = final_data
SCN5A['carrier_status'] = np.where(SCN5A['gene'] != 'SCN5A', 0, 1)
SCN5A = SCN5A.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
SCN5A['participant_type'] = np.where((SCN5A['participant_type']=='Rare Disease Proband'), 0, 1)
SCN5A['participant_type'] = np.where((SCN5A['participant_type']=='Cancer'), 2, SCN5A['participant_type'])
SCN5A['participant_phenotypic_sex'] = np.where((SCN5A['participant_phenotypic_sex']=='Male'), 0, 1)
SCN5A = SCN5A.dropna(how='all')
X = SCN5A.iloc[:,4:]
y = SCN5A[['carrier_status']]
 # building the model and fitting the data
log_reg_SCN5A = sm.Logit(y, X).fit(method='bfgs')
results_SCN5A = results_summary_to_dataframe(log_reg_SCN5A)
results_SCN5A = results_SCN5A.loc[results_SCN5A['pvals'] <= 0.05]
results_SCN5A.head(20)

# List with each logistic regression result dataframe per gene
data_frames = [results_MYBPC3, results_MYH7, results_TNNT2, results_TNNI3, 
               results_TPM1, results_MYL3, results_ACTC1, results_PRKAG2, 
               results_GLA, results_MYL2, results_LMNA, results_PKP2, 
               results_DSP, results_DSC2, results_TMEM43, results_DSG2, 
               results_KCNQ1, results_KCNH2, results_SCN5A]

# Dropping p-values, 95% CI columns
for df in data_frames:
    df.drop(['pvals', 'conf_lower', 'conf_higher'], axis=1, inplace=True)

# Genes of interest
genes = ['MYBPC3', 'MYH7', 'TNNT2', 'TNNI3', 'TPM1', 'MYL3', 'ACTC1', 'PRKAG2', 'GLA', 'MYL2', 'LMNA', 
         'PKP2', 'DSP', 'DSC2', 'TMEM43', 'DSG2', 'KCNQ1', 'KCNH2', 'SCN5A']

# Creating dataframe with all significant ORs per gene
results = pd.concat(data_frames, axis=1, levels = genes)

# Columns set to each gene
results.columns = genes

# Cleaning dataframe
results = results.reset_index()
results = results.rename(columns = {'index' : 'comorbidity'})

# Writing to file
results.to_csv("gene_comorbidities_logistic_reg.txt", index = False)

# Editing data to indicators per disease:
## carrier status indicating if participant is carrying particular disease (0= np, 1= yes)
## participant type: proband = 0, relative = 1, cancer = 2
## gender: 0 = male, 1 = female
## Subsetting into outcome (carrier status) and all covariates (age, gender, participant type, and all comorbidities)

# Building logistic regression and fitting to data
## Subsetting result to those with p-values of 0.05 or lower
# ## Per disease

HCM = final_data
HCM['carrier_status'] = np.where(((HCM['gene'] == 'MYBPC3') | 
                                  (HCM['gene'] == 'MYH7') | 
                                  (HCM['gene'] == 'TNNI3') | 
                                  (HCM['gene'] == 'TPM1') | 
                                  (HCM['gene'] == 'MYL3') | 
                                  (HCM['gene'] == 'ACTC1') | 
                                  (HCM['gene'] == 'MYL2') | 
                                  (HCM['gene'] == 'PRKAG2')), 1, 0)
HCM = HCM.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
HCM['participant_type'] = np.where((HCM['participant_type']=='Rare Disease Proband'), 0, 1)
HCM['participant_type'] = np.where((HCM['participant_type']=='Cancer'), 2, HCM['participant_type'])
HCM['participant_phenotypic_sex'] = np.where((HCM['participant_phenotypic_sex']=='Male'), 0, 1)
HCM = HCM.dropna(how='all')
X = HCM.iloc[:,4:]
y = HCM[['carrier_status']]
 # building the model and fitting the data
log_reg_HCM = sm.Logit(y, X).fit(method='bfgs')
results_HCM = results_summary_to_dataframe(log_reg_HCM)
results_HCM = results_HCM.loc[results_HCM['pvals'] <= 0.05]
results_HCM.head(20)

DCM = final_data
DCM['carrier_status'] = np.where(((DCM['gene'] == 'TNNT2') | 
                                  (DCM['gene'] == 'LMNA')), 1, 0)
DCM = DCM.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
DCM['participant_type'] = np.where((DCM['participant_type']=='Rare Disease Proband'), 0, 1)
DCM['participant_type'] = np.where((DCM['participant_type']=='Cancer'), 2, DCM['participant_type'])
DCM['participant_phenotypic_sex'] = np.where((DCM['participant_phenotypic_sex']=='Male'), 0, 1)
DCM = DCM.dropna(how='all')
X = DCM.iloc[:,4:]
y = DCM[['carrier_status']]
 # building the model and fitting the data
log_reg_DCM = sm.Logit(y, X).fit(method='bfgs')
results_DCM = results_summary_to_dataframe(log_reg_DCM)
results_DCM = results_DCM.loc[results_DCM['pvals'] <= 0.05]
results_DCM.head(20)

FD = final_data
FD['carrier_status'] = np.where((FD['gene'] == 'GLA'), 1, 0)
FD = FD.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
FD['participant_type'] = np.where((FD['participant_type']=='Rare Disease Proband'), 0, 1)
FD['participant_type'] = np.where((FD['participant_type']=='Cancer'), 2, FD['participant_type'])
FD['participant_phenotypic_sex'] = np.where((FD['participant_phenotypic_sex']=='Male'), 0, 1)
FD = FD.dropna(how='all')
X = FD.iloc[:,4:]
y = FD[['carrier_status']]
 # building the model and fitting the data
log_reg_FD = sm.Logit(y, X).fit(method='bfgs')
results_FD = results_summary_to_dataframe(log_reg_FD)
results_FD = results_FD.loc[results_FD['pvals'] <= 0.05]

ARVC = final_data
ARVC['carrier_status'] = np.where(((ARVC['gene'] == 'PKP2') |
                                   (ARVC['gene'] == 'DSP') |
                                   (ARVC['gene'] == 'DSC2') |
                                   (ARVC['gene'] == 'TMEM43') |
                                   (ARVC['gene'] == 'DSG2')), 1, 0)
ARVC = ARVC.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
ARVC['participant_type'] = np.where((ARVC['participant_type']=='Rare Disease Proband'), 0, 1)
ARVC['participant_type'] = np.where((ARVC['participant_type']=='Cancer'), 2, ARVC['participant_type'])
ARVC['participant_phenotypic_sex'] = np.where((ARVC['participant_phenotypic_sex']=='Male'), 0, 1)
ARVC = ARVC.dropna(how='all')
X = ARVC.iloc[:,4:]
y = ARVC[['carrier_status']]
 # building the model and fitting the data
log_reg_ARVC = sm.Logit(y, X).fit(method='bfgs')
results_ARVC = results_summary_to_dataframe(log_reg_ARVC)
results_ARVC = results_ARVC.loc[results_ARVC['pvals'] <= 0.05]

LQTS = final_data
LQTS['carrier_status'] = np.where(((LQTS['gene'] == 'KCNQ1') |
                                   (LQTS['gene'] == 'KCNH2') |
                                   (LQTS['gene'] == 'SCN5A')), 1, 0)
LQTS = LQTS.drop(columns=['participant_id','gene', 'participant_ethnic_category_simple'])
LQTS['participant_type'] = np.where((LQTS['participant_type']=='Rare Disease Proband'), 0, 1)
LQTS['participant_type'] = np.where((LQTS['participant_type']=='Cancer'), 2, LQTS['participant_type'])
LQTS['participant_phenotypic_sex'] = np.where((LQTS['participant_phenotypic_sex']=='Male'), 0, 1)
LQTS = LQTS.dropna(how='all')
X = LQTS.iloc[:,4:]
y = LQTS[['carrier_status']]
 # building the model and fitting the data
log_reg_LQTS = sm.Logit(y, X).fit(method='bfgs')
results_LQTS = results_summary_to_dataframe(log_reg_LQTS)
results_LQTS = results_LQTS.loc[results_LQTS['pvals'] <= 0.05]

# List of diseases
diseases = ['HCM', 'DCM', 'FD', 'ARVC', 'LQTS']

# # List with disease logistic regression result dataframes
data_frames = [results_HCM, results_DCM, results_FD, results_ARVC, results_LQTS]

# Dropping p-values, 95% CI columns
for df in data_frames:
    df.drop(['pvals', 'conf_lower', 'conf_higher'], axis=1, inplace=True)

# Merging disease dataframes and cleaning    
results = pd.concat(data_frames, axis=1, levels = genes)
results.columns = diseases
results = results.reset_index()
results = results.rename(columns = {'index' : 'comorbidity'})

# Writing to file
results.to_csv("disease_comorbidities_logistic_reg.txt", index = False)
