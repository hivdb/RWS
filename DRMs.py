import mysql.connector
import csv
import re
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


def combine_drms_from_two_columns(column1, column2):
    if (column1 == "" and column2 == ""):
        return ""
    elif column2 == "":
        return column1
    elif column1 == "":
        return column2
    else:
        return column1 + ", " + column2 


def combine_sort_drms_from_two_columns(column1, column2):
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


def sort_mutlist_bypos(mutlist_string):
    mutlist = mutlist_string.split(", ")
    sorted_mutations = sorted(mutlist, key=extract_position)
    sorted_mutlist_string = (", ").join(sorted_mutations)
    return sorted_mutlist_string


def extract_position(mutation):
    match = re.search(r'\d+', mutation)
    return int(match.group()) if match else 0


def get_scores(drug_class, mut_column):
    scored_muts = set()
    Score_file = f"{drug_class}_DataSets/{drug_class}_Scores.csv"
    with open(Score_file, mode='r', encoding='utf-8-sig') as file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            scored_muts.add(row[mut_column])
    return scored_muts


# sample_drms is a string. It may have more than one aa following the pos
# included_drms is a look-up list used to filter the sample_drms. It has just one aa following the pos
# Because M184VI is common, we want to prevent it from being a mixture so it is converted to M184V 
# Can this function be shortened in order to use sorted() with this as a key function
def filter_drms(sample_drms, included_drms):
    sample_drm_list = sample_drms.split(", ")
    filtered_list = []
    for sample_drm in sample_drm_list:
        if sample_drm == "M184VI":
            sample_drm = "M184V"
        pattern1 = r"\D(\d{1,3})(.+)"
        match1 = re.search(pattern1,sample_drm)
        if match1:
            sample_pos = match1.group(1)
            sample_mutations = match1.group(2)
        for included_drm in included_drms:
            pattern2 = r"\D(\d{1,3})(\D)"
            match2 = re.search(pattern2,included_drm)
            if match2:
                included_pos = match2.group(1)
                included_aa = match2.group(2)
            if sample_pos == included_pos:
                if included_aa in sample_mutations:
                    filtered_list.append(sample_drm) 
    return ", ".join(filtered_list)


# Insertions contain "_" and are not considered to be mixtures
def count_mixtures(drms):
    count = 0
    for drm in drms.split(', '):
        pattern = r"\D(\d{1,3})(.+)"
        match = re.search(pattern, drm)
        if match:
            aas = match.group(2)
            if "_" in aas:
                continue
            if len(aas) >=2:
                count += 1
    return count


# When a mutation contains a mixture of the consensus aa and a different aa, remove the consensus aa
# Otherwise the mutation is returned unchanged
def simplify_mutation(mutation):
    pattern = r'([A-Za-z])(\d+)([A-Za-z_]+)'
    match = re.match(pattern, mutation)
    if match:
        consensus, position, amino_acids = match.groups()
        if "_" in amino_acids: 
            return mutation
        elif consensus in amino_acids:
            simplified_amino_acids = amino_acids.replace(consensus, '')
            simplified_mutation = consensus + position + simplified_amino_acids
            return simplified_mutation
        else:
            return mutation


# This can be modified to provide multiple mutation when there is more than one 
#        non-consensus aa
# This function (1) counts the rows in the genotype-treatment files for a drug
# (2) counts the rows containing a DRM
# (3) creates a dictionary in which the key is a mutation and the value 
#            is the count of that mutation        
def create_drm_count_dict(file_path, column_name, scored_muts):
    mut_list = []
    num_isolates = 0
    with open(file_path, mode='r', encoding='utf-8') as file:
        csv_reader = csv.DictReader(file)
        header = next(csv_reader)
        for row in csv_reader:
            mut_list.append(row[column_name])
            num_isolates +=1 
    drm_counts = {}
    num_isolates_wdrm = 0
    for i in range(len(mut_list)):
        # Create a list of mutations
        sample_mut_list = mut_list[i].split(', ')
        sample_has_drm_flag = 0
        for j in range(len(sample_mut_list)):
            mutation = sample_mut_list[j]
            mutation = simplify_mutation(mutation)
            if mutation in scored_muts:
                sample_has_drm_flag = 1
                if mutation in drm_counts:
                    drm_counts[mutation] +=1
                else:
                    drm_counts[mutation]=1
        if sample_has_drm_flag == 1:
            num_isolates_wdrm +=1    
    return num_isolates, num_isolates_wdrm, drm_counts


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


def plot_drm_freqs(drug, num_isolates, num_isolates_wdrm, drm_freqs):
    #print("Sample DRMs:", drm_freq)
    #print(drms_to_show)
    # drm_freq_wout_cons = {}
    # for mut, freq in drm_freq.items(): 
    #     new_mut = mut[1:]
    #     drm_freq_wout_cons[new_mut] = freq

    # create a new dictionary with only muts_to_show
    # drm_freq_to_show = {element: 0 for element in drms_to_show}
    # for drm in drm_freq:
    #     if drm in drms_to_show: 
    #         drm_freq_to_show[drm] = drm_freq_to_show[drm]
    #print(drm_freq_to_show)    
    #drms_sorted_bypos = dict(sorted(drm_freq_wout_cons.items(), key=lambda item: item[0]))
    #drms_sorted_bypos = dict(sorted(drm_freq_to_show.items(), key=lambda item: sort_mut_bypos(item[1])))
    #drms_sorted_bypos = drm_freq_to_show

    mutations = list(drm_freqs.keys())
    pcnts = list(drm_freqs.values())
    plt.figure(figsize=(10, 3))
    plt.title(f"{drug} NumIsolates:{num_isolates} NumIsolatesWDRM: {num_isolates_wdrm} ")
    plt.bar(mutations, pcnts)
    plt.xlabel('DRMs')
    plt.xticks(rotation=90)
    plt.subplots_adjust(bottom=0.3)
    plt.ylabel('Percent')
    plt.ylim(0,70)
    plt.grid(axis='y')
    plt.savefig(drug+".png")


def is_iterable(obj):
    try:
        iter(obj)
        return True
    except TypeError:
        return False


def remove_cons_from_mutlist(mut_list):
    drms = mut_list.split(", ")
    drms_wout_cons = []
    for drm in drms:
        drms_wout_cons.append(drm[1:])
    return ', '.join(drms_wout_cons)


def plot_phenotypes(drug, phenotype_data, min_count=2, max_drms=7, \
                    min_median=3.0, xlim = 0.5):
    df = pd.DataFrame(phenotype_data, columns = ['Pattern', 'Count', 'NumDRMs', 'Folds', 'Median', \
                                                  'Min', 'Max', 'Q1', 'Q3'])
    df = df[df['Pattern'] != ""]
    df['Pattern'] = df['Pattern'].apply(remove_cons_from_mutlist)
    df = df[df['Count'] >= min_count]
    df = df[df['Median'] >= min_median]
    df = df[df['NumDRMs'] <= max_drms]
    df['Rank'] =df['Median'].rank(ascending=False, method='first')
    num_rows = len(df)
    #print(df)
    rows_to_add = []
    for index, row in df.iterrows():
        rank = row['Rank']
        pattern = row['Pattern']
        string_list_folds = row['Folds']
        median = row['Median']
        q1 = row["Q1"]
        q3 = row["Q3"]
        actual_list_folds = ast.literal_eval(string_list_folds)
        if not is_iterable(actual_list_folds):
            actual_list_folds = [actual_list_folds]
        for fold in actual_list_folds:
            if fold <xlim:
                fold = xlim
            if fold >128:
                fold=128
            new_row = (rank, pattern, fold)
            rows_to_add.append(new_row)

    df_for_plot = pd.DataFrame(rows_to_add)
    print(df_for_plot)
    column_names = ('Rank', 'Pattern', 'Fold')
    df_for_plot.columns = column_names
    plt.figure(figsize=(10, num_rows/6))
    plt.subplots_adjust(top=0.95)
    plt.subplots_adjust(left=0.3)
    sns.scatterplot(data=df_for_plot, x='Fold', y='Rank', color="white")
    x_tick_positions = [1, 2, 4, 8, 16, 32, 64, 128]
    x_tick_labels = ["1", "2", "4", "8", "16", "32", "64", "128"]
    ylim = plt.ylim()
    plt.ylim(ylim[::-1])
    plt.yticks(fontsize=8)
    plt.ylabel('')
    plt.yticks(ticks=df_for_plot['Rank'], labels=df_for_plot['Pattern'])
    plt.grid(axis='y')
    plt.xscale('log', base=2)
    plt.xlim(left=0.6)
    plt.xticks(ticks=x_tick_positions, labels=x_tick_labels)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=12)
    #plt.show()
    #plt.savefig(drug+".phenotypes.png")

    df_medians = df[['Pattern', 'Count', 'Rank', 'Median']] 
    df_q1 = df[['Pattern', 'Count', 'Rank', 'Q1']] 
    df_q3 = df[['Pattern', 'Count', 'Rank', 'Q3']] 
    df_medians.loc[df_medians["Median"] > 128, "Median"] = 128
    plt.subplots_adjust(top=0.95)
    title = f"{drug} susc (Min #results/pattern:{min_count} \n" + \
            f" Lowest median:{min_median} Max DRMs:{max_drms})"
    plt.title(title, fontsize=16, fontweight='bold') 
    plt.subplots_adjust(left=0.3)
    sns.scatterplot(data=df_medians, x='Median', y='Rank', color="blue", label="Median")
    sns.scatterplot(data=df_q1, x='Q1', y='Rank', color="blue", marker="<", label="IQ1" )
    sns.scatterplot(data=df_q3, x='Q3', y='Rank', color='blue', marker=">", label="IQ3")
    plt.ylabel('')
    plt.xlabel('Fold Reduced Susceptibility', fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=12)
    #plt.show()
    plt.savefig(drug+".phenotypes_medians.png")
