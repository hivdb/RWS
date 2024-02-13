import csv
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

drms_to_show = {"NRTI": ["M41L", "A62V", "K65R", "K65N", "D67N", "T69_", "K70R", "K70E", "K70Q", "L74V", "L74I",
                        "Y115F", "Q151M", "M184V", "M184I", "L210W", "T215Y", "T215F", "K219E", "K219Q"],
                "NNRTI": ["A98G", "L100I", "K101E", "K101P", "K103N", "V106A", "V106M", "E138K", "V179D", 
                          "Y181C", "Y181I", "Y181V", "Y188C", "Y188L", "G190A", "G190S", "G190E", "H221Y", 
                          "P225H", "F227L", "F227C", "M230L", "Y318F"],
                "INSTI": ["H51Y", "T66I", "T66K", "E92Q", "T97A", "G118R", "F121Y", "E138K", "E138A", 
                          "G140S", "G140A", "S147G", "Q148H", "Q148R", "Q148K", "S153Y", "N155H", "R263K"],
                "PI": ["L10F", "L24I", "V32I", "L33F", "M46I", "M46L", "I47V", "I47A", "G48V", 
                       "I50L", "I50V", "F53L", "I54V", "I54L", "I54M", "G73S", "T74P", "L76V", "V82A", 
                       "V82T", "V82F", "I84V", "N88S", "L89V", "L89T", "L90M"]}

from DRMs import get_scores, create_drm_count_dict, plot_drm_freqs

csv_files = {
   # "NRTI": ["ABC_3TC", "AZT_3TC", "TDF_XTC"], 
   # "NNRTI": ["EFV", "RPV", "DOR"],
    "INSTI": ["CAB", "DTG", "EVG", "RAL"]}
   # "PI": ["ATV", "LPV"]}

for drug_class, files in csv_files.items():
    dir = drug_class + "_DataSets"
    if drug_class in ["NRTI", "NNRTI"]: 
        col_name = drug_class + "DRMs" 
    elif drug_class == "INSTI": 
        col_name = "INIMajorDRMs"
    else:
        col_name = "CompMutList" 

    scored_muts = get_scores(drug_class, "Rule")
    drms_to_show = drms_to_show[drug_class]
    
    for drug in files: 
        path = dir + "/" + drug + ".csv"
        
        # drm_counts is a dictionary in which the keys are scored drms and the values
        # are the number of occurrences of the drm
        (num_isolates, num_isolates_wdrm, drm_sample_counts) = create_drm_count_dict(path, col_name, 
                                                                              scored_muts)
             
        # create a new dictionary in which drms_to_show which are not in drm_counts are added
        # and assigned a count of 0
        drm_all_counts = {k: drm_sample_counts.get(k, 0) for k in 
                          set(drms_to_show + list(drm_sample_counts.keys()))}
        df = pd.DataFrame(list(drm_all_counts.items()), columns =['DRM', 'Count'])
        
        # Filter the df so that is only has rows with DRMs that are in drms_to_show
        df = df[df['DRM'].isin(drms_to_show)]
        df['Pcnt'] = (df['Count'] / num_isolates_wdrm) * 100
        df['Pcnt'] = df['Pcnt'].round(1)
        df[['Pos', 'Mut']] = df['DRM'].str.extract(r'(\d+)([A-Za-z])')
        df['Pos'] = pd.to_numeric(df['Pos'])
        df = df.sort_values(by=['Pos', 'Mut'], ascending=[True, True])
        plot_drm_freqs(drug, num_isolates, num_isolates_wdrm, df)
