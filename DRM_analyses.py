import csv
import re
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
        (num_isolates, num_isolates_wdrm, drm_counts) = create_drm_count_dict(path, col_name, 
                                                                              scored_muts)

        sorted_drm_freq = dict(sorted(drm_counts.items(), key=lambda x: x[1], reverse=True))

        for drm, count in sorted_drm_freq.items():
            pcnt = (count/num_isolates_wdrm) * 100
            rounded_pcnt = round(pcnt, 1)
            sorted_drm_freq[drm] = rounded_pcnt

        sorted_drm_freq_to_show ={drm: sorted_drm_freq[drm] for drm in drms_to_show if drm in sorted_drm_freq}

        plot_drm_freqs(drug, num_isolates, num_isolates_wdrm, sorted_drm_freq_to_show)
