import csv
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import statistics
import numpy as np

drug_class_drugs = {"NRTI": ["ABC", "TDF"], 
                    "NNRTI": ["EFV"],
                    "PI": ["ATV", "DRV", "LPV"],
                    "INSTI": ["BIC", "CAB", "DTG"]}

drms_to_show = {"NRTI": ["M41L", "A62V", "K65R", "K65N", "D67N", "T69_", "K70R", "K70E", "K70Q", "L74V", "L74I",
                        "Y115F", "A151M", "M184V", "M184I", "L210W", "T215Y", "T215F", "K219E", "K219Q"],
                "NNRTI": ["A98G", "L100I", "K101E", "K101P", "K103N", "K103S", "V106A", "V106M", "E138A", 
                          "E138K", "V179D", "Y181C", "Y181I", "Y181V", "Y188C", "Y188L", "G190A", "G190S", 
                          "G190E", "H221Y", "P225H", "F227L", "F227C", "M230L", "Y318F"],
                "INSTI": ["H51Y", "T66I", "T66K", "E92Q", "T97A", "G118R", "F121Y", "E138K", "E138A", 
                          "G140S", "G140A", "S147G", "Q148H", "Q148R", "Q148K", "S153Y", 
                          "N155H", "R263K"],
                "PI": ["L10F", "L24I", "V32I", "L33F", "M46I", "M46L", "I47V", "I47A", "G48V", "G48M", 
                       "I50L", "I50V", "F53L", "I54V", "I54L", "I54M", "I54S", "I54T", "I54A", "G73S", "T74P", 
                        "L76V", "V82A", "V82T", "V82F", "I84V", "N88S", "L89V", "L89T", "L90M"]}

from DRMs import (combine_drms_from_two_columns, sort_mutlist_bypos, filter_drms,
                  count_mixtures, plot_phenotypes)

# Main program phenotype file
# Columns: RefID, IsolateID, IsolateName, Species, Type (Clinical vs Lab), PtID, Method
# DOR, DORFoldMatch, EFV, EFVFoldMatch, ETR, ETRFoldMatch, NVP, NVPFoldMatch, RPV, RPVFoldMatch, 
# CompleteMutationListAvailable, NRTIDRMs, NNRTIDRMs, NonDRMs, P1 ... P300 
drug_class = "PI"
dir = drug_class + "_DataSets"
file = drug_class + "_Phenotypes.csv"
csv_file = dir + "/" + file
#print(csv_file)
pcnt_mix_disqualified = 0.3

column_to_read = []
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

for drug in drug_class_drugs[drug_class]:
    phenotypes = []
    drm_pattern_counts = {}
    drm_pattern_folds = {}
    drm_patterns = {}
    num_isolates = 0
    with open(csv_file, mode='r', encoding='utf-8-sig') as file:      
        reader = csv.DictReader(file)
        selected_data = [] 
        for row in reader:
            selected_row = {column: row[column] for column in columns_to_read}
            #print("selected_row", selected_row)
            selected_data.append(selected_row) 
            if drug_class == "PI":
                all_drms = combine_drms_from_two_columns(selected_row["PIMajorDRMs"], 
                                                         selected_row["PIMinorDRMs"])
            elif drug_class == "INSTI":
                all_drms = combine_drms_from_two_columns(selected_row["INIMajorDRMs"], 
                                                         selected_row["INIMinorDRMs"])
            else:
                all_drms = selected_row[drm_column]
            all_drms = sort_mutlist_bypos(all_drms)       
            if selected_row["Method"] != "PhenoSense" or selected_row[drug] == "" or all_drms == "":
                continue          
            filtered_drms = filter_drms(all_drms, drms_to_show[drug_class])         
            phenotype = float(selected_row[drug].replace(",", ""))
            num_drms = len(filtered_drms.split(', '))
            num_mixtures = count_mixtures(filtered_drms)            
            if (num_mixtures/num_drms) > pcnt_mix_disqualified:
                continue
            num_isolates +=1    
            phenotypes.append([filtered_drms, drug, phenotype])     
        #for item in phenotypes:
        #    print(item)
          
    for (filtered_drms, drug, fold) in phenotypes:
        if filtered_drms in drm_patterns:
            drm_patterns[filtered_drms]["count"] +=1
            drm_patterns[filtered_drms]["fold"].append(fold)
        else:
            drm_patterns[filtered_drms] = {"count": 1, "fold": [fold]}            

    drms_phenotype_table = []       
    for filtered_drms in drm_patterns.keys():
        count = drm_patterns[filtered_drms]["count"]
        num_drms = len(filtered_drms.split(', '))
        median_fold = statistics.median(drm_patterns[filtered_drms]["fold"])
        min_fold = min(drm_patterns[filtered_drms]["fold"])
        max_fold = max(drm_patterns[filtered_drms]["fold"])
        q1 = np.percentile(drm_patterns[filtered_drms]["fold"],25)
        q3 = np.percentile(drm_patterns[filtered_drms]["fold"],75)
        string_folds = str(drm_patterns[filtered_drms]["fold"])
        string_folds = string_folds[1:-1]
        #print(type(string_folds), "  ", string_folds)
        drms_data = (filtered_drms, count, num_drms, string_folds, median_fold, min_fold, max_fold, q1, q3)
        drms_phenotype_table.append(drms_data)
    
    num_patterns = len(drms_phenotype_table)
    sorted_drms_table_by_count = sorted(drms_phenotype_table, key=lambda x:x[1], reverse=True)
    sorted_drms_table_by_fold = sorted(drms_phenotype_table, key=lambda x:x[4], reverse=True)
    #for item in sorted_drms_table_by_fold:
    #    print(item)

    plot_phenotypes(drug, sorted_drms_table_by_fold)

    file_name = drug + "_drm_patterns.csv"
    with open(file_name, 'w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)
        for row in sorted_drms_table_by_fold:
            writer.writerow(row)
    