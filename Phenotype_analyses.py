import csv
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import statistics
import numpy as np

pd.set_option('display.max_rows', 300)

drug_class_drugs = {"NRTI": ["TDF"], 
                    "NNRTI": ["EFV"],
                    "PI": ["ATV"],
                    "INSTI": ["BIC", "CAB", "DTG"]}

drms_to_show = {"NRTI": ["M41L", "A62V", "K65R", "K65N", "D67N", "T69i", "K70R", "K70E", "K70Q", "L74V", "L74I",
                        "Y115F", "A151M", "M184V", "M184I", "L210W", "T215Y", "T215F", "K219E", "K219Q"],
                "NNRTI": ["A98G", "L100I", "K101E", "K101P", "K103N", "K103S", "V106A", "V106M", "E138A", 
                          "E138K", "V179D", "Y181C", "Y181I", "Y181V", "Y188C", "Y188L", "G190A", "G190S", 
                          "G190E", "H221Y", "P225H", "F227L", "F227C", "M230L", "Y318F"],
                "INSTI": ["H51Y", "T66I", "T66K", "E92Q", "T97A", "G118R", "F121Y", "E138K", "E138A", 
                          "G140S", "G140A", "S147G", "Q148H", "Q148R", "Q148K", "S153Y", 
                          "N155H", "R263K"],
                "PI": ["L10F", "L24I", "V32I", "L33F", "M46I", "M46L", "I47V", "I47A", "G48V", "G48M", 
                       "I50L", "I50V", "F53L", "I54V", "I54L", "I54M", "I54S", "I54T", "I54A", "G73S", "T74P", 
                        "L76V", "V82A", "V82T", "V82F", "I84V", "N88S", "N88D", "L89V", "L89T", "L90M"]}

from DRMs import (combine_drms_from_two_columns, format_indels, sort_mutlist_bypos, filter_drms,
                  count_mixtures, plot_phenotypes)

# Main program phenotype file
# Columns: RefID, IsolateID, IsolateName, Species, Type (Clinical vs Lab), PtID, Method
# DOR, DORFoldMatch, EFV, EFVFoldMatch, ETR, ETRFoldMatch, NVP, NVPFoldMatch, RPV, RPVFoldMatch, 
# CompleteMutationListAvailable, NRTIDRMs, NNRTIDRMs, NonDRMs, P1 ... P300 
drug_class = "NRTI"
dir = drug_class + "_DataSets"
file = drug_class + "_Phenotypes.csv"
csv_file = dir + "/" + file
pcnt_mix_disqualified = 0.3
df= pd.read_csv(csv_file, encoding='utf-8-sig')

if drug_class == "NNRTI":
    columns_to_read = ["Method", "NNRTIDRMs", "DOR", "EFV", "RPV"]
    drm_column = "NNRTIDRMs"
elif drug_class == "NRTI":
    columns_to_read = ["Method", "NRTIDRMs", "ABC", "AZT", "TDF"]
    drm_column = "NRTIDRMs"
elif drug_class == "PI":
    columns_to_read = ["Method", "PIMajorDRMs", "PIMinorDRMs", "ATV", "DRV", "LPV"]
    drm_column = ""
elif drug_class == "INSTI":
    columns_to_read = ["Method", "INIMajorDRMs", "INIMinorDRMs", "BIC", "CAB", "DTG"]
    drm_column = ""

# For each drug, place all of the DRMs for its drug class into 'AllDRMs'  
for drug in drug_class_drugs[drug_class]:
    df = df[columns_to_read].copy()
    if drug_class == "PI":
        df['AllDRMs'] = df.apply(lambda row: combine_drms_from_two_columns(row['PIMajorDRMs'], row['PIMinorDRMs']), axis=1)
    elif drug_class == "INSTI":
        df['AllDRMs'] = df.apply(lambda row: combine_drms_from_two_columns(row['INIMajorDRMs'], row['INIMinorDRMs']), axis=1)
    else:
        df['AllDRMs'] = df[drm_column]
       
    # Filter df to contain only rows with DRMs and drug susc results performed by PhenoSense assy  
    df = df.loc[(df['Method'] == "PhenoSense") & (df[drug].notna()) & (df['AllDRMs'] != ""), :].copy()
    df = df.loc[(df['AllDRMs'].notna()), :].copy()
    
    #Convert insertions such as "T69S_SS" to T69i and Deletions such as T69~ to T69d
    df['AllDRMs'] = df['AllDRMs'].astype(str)
    df['AllDRMs'] = df['AllDRMs'].apply(format_indels)
    #print(df.head(100))

    # Add a column that shows only those DRMs that are in drms_to_show
    df.loc[:, 'DRMsToShow'] = df['AllDRMs'].apply(lambda drms: filter_drms(drms, drms_to_show[drug_class]))
    df = df.loc[(df['DRMsToShow'] != ""), :].copy()
     
    if df[drug].dtype == 'object':
        df.loc[:, 'Fold'] = df[drug].str.replace(",", "").astype(float)
    else:
        df.loc[:, 'Fold'] = df[drug].astype(float)

    # Further filtering based on mixture percentage
    df.loc[:, 'NumDRMs'] = df['DRMsToShow'].apply(lambda x: len(x.split(', ')))
    df.loc[:, 'NumMixtures'] = df['DRMsToShow'].apply(count_mixtures)
    df['PcntMix'] = df['NumMixtures']/df['NumDRMs']
    df = df[df['PcntMix'] < 0.3]

    # Group by 'filtered_drms' and aggregate
    drm_group = df.groupby('DRMsToShow').agg(
        Count=('DRMsToShow', 'size'),  # Count the occurrences of each DRM pattern
        Folds=('Fold', list)  # Collect all phenotype values in a list for each DRM pattern
    ).reset_index()
    #print(drm_group.head(300))

    # Convert the aggregated DataFrame to the desired dictionary format
    drm_patterns = drm_group.set_index('DRMsToShow').T.to_dict('list')
    #print("\n\n", drm_patterns)

   # Convert drm_patterns dictionary to DataFrame
    drm_patterns_df = pd.DataFrame.from_dict(drm_patterns, orient='index', 
                                    columns=['Count', 'Folds']).reset_index()
    drm_patterns_df.rename(columns={'index': 'DRMsToShow'}, inplace=True)
    #print("\n\n", drm_patterns_df)

    # Calculate the necessary statistics
    drm_patterns_df['Median'] = drm_patterns_df['Folds'].apply(np.median)
    drm_patterns_df['Min'] = drm_patterns_df['Folds'].apply(min)
    drm_patterns_df['Max'] = drm_patterns_df['Folds'].apply(max)
    drm_patterns_df['Q1'] = drm_patterns_df['Folds'].apply(lambda x: np.percentile(x, 25))
    drm_patterns_df['Q3'] = drm_patterns_df['Folds'].apply(lambda x: np.percentile(x, 75))
    drm_patterns_df['StringFolds'] = drm_patterns_df['Folds'].apply(lambda x: str(x)[1:-1])
    drm_patterns_df['NumDRMs'] = drm_patterns_df['DRMsToShow'].apply(lambda x: len(x.split(', ')))
    drm_patterns_df = drm_patterns_df[['DRMsToShow', 'Count', 'NumDRMs', 'StringFolds', 'Median',
                                       'Min', 'Max', 'Q1', 'Q3']]
    
    # Sort by 'median_fold' in descending order
    drm_patterns_df.sort_values(by='Median', ascending=False, inplace=True)
    #print("\n\n", drm_patterns_df)
    
    # Convert DataFrame to the list of tuples if needed for the plot_phenotypes function
    drm_patterns_table_by_fold = [tuple(x) for x in drm_patterns_df.to_numpy()]
    
    plot_phenotypes(drug, drm_patterns_table_by_fold)

    file_name = drug + "_drm_patterns.csv"
    #with open(file_name, 'w', newline='', encoding='utf-8') as file:
    #    writer = csv.writer(file)
    #    for row in drm_patterns_table_by_fold:
    #        writer.writerow(row)
    