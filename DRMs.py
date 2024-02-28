import mysql.connector
import csv
import re
import logging
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import ast

# Configure database connection
# Need to use Philip's code to set up the environment at the command line: bash phpmyadmin
config = {
    'user': 'rshafer',
    'password': 'rshafer',
    'host': '10.77.6.244',
    'database': 'HIVDB2',
    'raise_on_warnings': True
}

tams = {41: ["*"], 67: ["*"], 70: ["R"], 210: ["*"], 215: ["*"], 219: ["*"]}

with open('app.log', 'a') as log_file:
    print("This is an info message", file=log_file)
logging.basicConfig(filename='app.log', level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')


def is_iterable(obj):
    try:
        iter(obj)
        return True
    except TypeError:
        return False


def is_tam(drm):
    pattern = r'([A-Za-z])(\d+)([A-Za-z]+)'
    match = re.match(pattern, drm)
    if match:
        cons, pos, mut = match.groups()
        pos = int(pos)
        if pos not in tams: 
            return False
        if pos == 70 and mut != "R":
            return False
        return True
    return False

#combines mutations in two strings from  a df into one string
#and sorts them by their positions
def combine_drms_from_two_columns(column1, column2):
    # Handle None or NaN values by converting them to empty strings
    column1 = "" if pd.isnull(column1) else column1
    column2 = "" if pd.isnull(column2) else column2
    new_string_list = ""
    if (column1 == "" and column2 == ""):
        new_string_list = ""
    elif column2 == "":
        new_string_list = column1
    elif column1 == "":
        new_string_list = column2
    else:
        new_string_list = column1 + ", " + column2 
        new_string_list = sort_mutlist_bypos(new_string_list)
    return new_string_list


# Accepts a string with a list of comma-separated mutattion
# returns a string in which the mutations are sorted by position        
def sort_mutlist_bypos(mutlist_string):
    mutlist = mutlist_string.split(", ")
    sorted_mutations = sorted(mutlist, key=extract_position)
    sorted_mutlist_string = (", ").join(sorted_mutations)
    return sorted_mutlist_string


def extract_position(mutation):
    match = re.search(r'\d+', mutation)
    return int(match.group()) if match else 0


# sample_drms is a string. It may have more than one aa following the pos
# included_drms is a look-up list used to filter the sample_drms. It has just one aa following the pos
# Because M184VI is common, we want to prevent it from being a mixture so it is converted to M184V 
# Can this function be shortened in order to use sorted() with this as a key function
# HOW CAN THIS BE SIMPLIFIED???
def filter_drms(sample_drms, drms_to_show):
    #print ("\n\nIn filter_drms, sample_drms", sample_drms)
    sample_drm_list = sample_drms.split(", ")
    filtered_list = []
    for sample_drm in sample_drm_list:
        if sample_drm == "M184VI":
            sample_drm = "M184V"
        pattern1 = r"(\D)(\d{1,3})(.+)"
        match1 = re.search(pattern1,sample_drm)
        if match1:
            cons = match1.group(1)
            sample_pos = match1.group(2)
            sample_mutations = match1.group(3)
            if "_" in sample_mutations:
                sample_drm = cons + sample_pos + "i"
                sample_drms = f"T69i, {sample_drms}"
                print(sample_drms)
        for drm in drms_to_show:
            pattern2 = r"\D(\d{1,3})(\D)"
            match2 = re.search(pattern2, drm)
            if match2:
                included_pos = match2.group(1)
                included_aa = match2.group(2)
            if sample_pos == included_pos:
                if int(sample_pos) == 69:
                    print("Sample DRM:" , sample_drm, "Included aa: ", included_aa)
                if included_aa in sample_mutations:

                    filtered_list.append(sample_drm) 
                if "i" in sample_drm:
                    print(sample_drm)
    #if "T69i" in filtered_list:
    #    print(filtered_list)
    return ", ".join(filtered_list)
    new_list = ", ".joint(filtered_list)
    if "T69i" in new_list:
        print("Filtered List: ", filtered_list)


# Insertions contain "_" and are not considered to be mixtures
# Consider indicating insertions with "Ins"
def count_mixtures(drms):
    count = 0
    if "T69i" in drms:
        print(drms)
    for drm in drms.split(', '):
        pattern = r"\D(\d{1,3})(.+)"
        match = re.search(pattern, drm)
        if match:
            aas = match.group(2)
            if aas == "i" in aas:
                print ("In count_mixtures", drm)
            if len(aas) >=2:
                count += 1
    return count


# When a mutation contains a mixture of the consensus aa and a different aa, remove the consensus aa
# Otherwise the mutation is returned unchanged
def simplify_mutation(mutation):
    pattern = r'([A-Za-z])(\d+)([A-Za-z_~]+)'
    match = re.match(pattern, mutation)
    if match:
        consensus, position, amino_acids = match.groups()
        if "_" in amino_acids: 
            new_mutation = consensus + position + "i"
            return new_mutation
        elif amino_acids == "~":
            new_mutation = consensus + position + "d"
            return new_mutation
        elif consensus in amino_acids:
            new_amino_acids = amino_acids.replace(consensus, '')
            new_mutation = consensus + position + new_amino_acids
            return new_mutation
        else:
            return mutation


def create_drm_count_dict(list_of_sample_muts, scored_muts_list):
    num_isolates_wdrm = 0
    drm_counts = {}
    for i in range(len(list_of_sample_muts)):
          sample_muts = list_of_sample_muts[i].split(', ')
          sample_has_drm_flag = 0
          for j in range(len(sample_muts)):
              mutation = sample_muts[j]
              mutation = simplify_mutation(mutation)
              if mutation in scored_muts_list:
                  sample_has_drm_flag = 1
                  if mutation in drm_counts:
                      drm_counts[mutation] +=1
                  else:
                      drm_counts[mutation]=1
          if sample_has_drm_flag == 1:
              num_isolates_wdrm +=1
    return num_isolates_wdrm, drm_counts
          

def create_drm_pattern_count_dict(file_path, column_name, scored_muts, muts_to_show):
    mut_list = []
    num_isolates = 0
    with open(file_path, mode='r', encoding='utf-8') as file:
        csv_reader = csv.DictReader(file)
        header = next(csv_reader)
        for row in csv_reader:
            mut_list.append(row[column_name])
            num_isolates +=1
    drm_patterns = {}
    for i in range(len(mut_list)):
        sample_mut_list = mut_list[i].split(', ')
        pattern = []
        for j in range(len(sample_mut_list)):
            mutation = sample_mut_list[j]
            mutation = simplify_mutation(mutation)
            if (mutation in scored_muts) and (mutation[1:] in muts_to_show):
                pattern.append(mutation)
        pattern = tuple(pattern)
        if pattern in drm_patterns:
            drm_patterns[pattern] +=1
        else:
            drm_patterns[pattern]=1
    sorted_drm_patterns = sorted(drm_patterns.items(), key=lambda x: x[1], reverse=True)
    sorted_drmpatterns_pcnt = [(mutation, round(100* freq/num_isolates, 1)) for mutation, freq in sorted_drm_patterns]
    return sorted_drmpatterns_pcnt


def sort_mut_bypos(mut):
    numeric_part = mut[:-1]
    return int(numeric_part)


def extract_position(mut):
    match = re.search(r'\d+', mut)
    return int(match.group()) if match else None


def plot_drm_freqs(drug, num_isolates, num_isolates_wdrm, drm_freqs_df):
    drms = list(drm_freqs_df["DRM"])
    pcnts = list(drm_freqs_df["Pcnt"])
    plt.figure(figsize=(10, 3))
    plt.title(f"{drug} NumIsolates:{num_isolates} NumIsolatesWDRM: {num_isolates_wdrm} ")
    plt.bar(drms, pcnts)
    plt.xlabel('DRMs')
    plt.xticks(rotation=90)
    plt.subplots_adjust(bottom=0.3)
    plt.ylabel('Percent')
    plt.ylim(0,60)
    plt.grid(axis='y')
    plt.savefig(drug+".png")


def remove_cons_from_mutlist(mut_list):
    drms = mut_list.split(", ")
    drms_wout_cons = []
    for drm in drms:
        drms_wout_cons.append(drm[1:])
    return ', '.join(drms_wout_cons)


def plot_phenotypes(drug, phenotype_data, min_count=1, max_drms=7, 
                    min_median=4.0, xlim = 0.5):
    df = pd.DataFrame(phenotype_data, columns = ['Pattern', 'Count', 'NumDRMs', 'Folds', 'Median', \
                                                  'Min', 'Max', 'Q1', 'Q3'])
    df['Pattern'] = df['Pattern'].apply(remove_cons_from_mutlist)
    df = df[df['Count'] >= min_count]
    df = df[df['Median'] >= min_median]
    df = df[df['NumDRMs'] <= max_drms]
    df['Rank'] =df['Median'].rank(ascending=True, method='first')
    num_rows = len(df)

    medians_df = df[['Pattern', 'Count', 'Rank', 'Median']] 
    q1_df = df[['Pattern', 'Count', 'Rank', 'Q1']] 
    q3_df = df[['Pattern', 'Count', 'Rank', 'Q3']] 
    medians_df.loc[medians_df["Median"] > 128, "Median"] = 128
    q1_df.loc[q1_df["Q1"] > 128, "Q1"] = 128
    q3_df.loc[q3_df["Q3"] > 128, "Q3"] = 128
    title = (
        f"{drug} susc (min #results/pattern:{min_count})"  
        #f"\nlowest median:{min_median} max #DRMs:{max_drms})"
        )
    x_tick_positions = [1, 2, 4, 8, 16, 32, 64, 128]
    x_tick_labels = ["1", "2", "4", "8", "16", "32", "64", "128"]
    fold_cutoffs = [5.2]
    plt.figure(figsize=(10, num_rows/6))
    plt.title(title, fontsize=16, fontweight='bold') 
    sns.scatterplot(data=medians_df, x='Median', y='Rank', color="blue", label="Median")
    sns.scatterplot(data=q1_df, x='Q1', y='Rank', color="blue", marker="<", label="IQ1" )
    sns.scatterplot(data=q3_df, x='Q3', y='Rank', color='blue', marker=">", label="IQ3")
    plt.subplots_adjust(top=0.95)
    plt.subplots_adjust(left=0.3)
    plt.ylabel('')
    plt.yticks(ticks=medians_df['Rank'], labels=medians_df['Pattern'])
    plt.yticks(fontsize=10)
    plt.grid(axis='y')
    plt.xscale('log', base=2)
    plt.xlabel('Fold Reduced Susceptibility', fontsize=14)
    plt.xticks(ticks=x_tick_positions, labels=x_tick_labels)
    plt.xticks(fontsize=14)
    plt.xlim(left=0.6)
    for fold in fold_cutoffs:
        plt.axvline(x=fold, color='k', linestyle=':')
    #plt.show()
    plt.savefig(drug+".phenotypes_medians.png")
