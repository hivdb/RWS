import mysql.connector
import csv
import re
import seaborn as sns
import matplotlib.pyplot as plt

# Configure database connection
# Need to use Philip's code to set up the environment at the command line: bash phpmyadmin
config = {
    'user': 'rshafer',
    'password': 'rshafer',
    'host': '10.77.6.244',
    'database': 'HIVDB2',
    'raise_on_warnings': True
}

muts_to_show = {"NRTI": ["41L", "62V", "65R", "65N", "67N", "70R", "70E", "70Q", "74V", "74I", \
                        "115F", "151M", "184V", "184I", "210W", "215Y", "215F", \
                        "219E", "219Q"],
                "NNRTI": ["98G", "100I", "101E", "101P", "103N", "106A", "106M", "138K", "179D", \
                        "181C", "188C", "188L", "190A", "190S", "221Y", "225H", "227L", "227C", "230L"],
                "INSTI": ["51Y", "66I", "66K", "92Q", "97A", "118R", "121Y", "138K", "138A", "140S", "140A", \
                          "143C", "143R", "147G", "148H", "148R", "148K", "153Y", "155H", "263K"],
                "PI": ["10F", "24I", "32I", "33F", "46I", "47V", "47A", "48V", "50L", "50V", \
                       "54V", "54L", "54M", "73S", "74P", "76V", "82A", "82T", "82F", "84V", "88S", \
                        "89V", "89T", "90M"]}

from DRMs import *

# Main program
csv_files = {
    #"NRTI": ["ABC_3TC", "AZT_3TC", "TDF_XTC"], 
    #"NNRTI": ["EFV", "RPV", "DOR"],
    #"INSTI": ["CAB", "DTG", "EVG", "RAL"],
    "PI": ["ATV", "LPV"]}

for drug_class, files in csv_files.items():
    dir = drug_class + "_DataSets"
    if drug_class in ["NRTI", "NNRTI"]: 
        col_name = drug_class + "DRMs" 
    elif drug_class == "INSTI": 
        col_name = "INIMajorDRMs"
    else:
        col_name = "CompMutList" 

    scored_muts = get_scores(drug_class, "Rule")
    
    
    for drug in files: 
        path = dir + "/" + drug + ".csv"
        (num_isolates, num_isolates_wdrm, drm_counts) = create_drm_count_dict(path, col_name, scored_muts)

        sorted_drm_freq = dict(sorted(drm_counts.items(), key=lambda x: x[1], reverse=True))
        
        for drm, count in sorted_drm_freq.items():
            pcnt = (count/num_isolates_wdrm) * 100
            rounded_pcnt = round(pcnt, 1)
            sorted_drm_freq[drm] = rounded_pcnt
        

        plot_drm_freqs(drug, num_isolates, num_isolates_wdrm, sorted_drm_freq, muts_to_show[drug_class])
          
        print("\n\n")