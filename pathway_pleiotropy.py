"""
pathway_pleiotropy.py
This script analyzes gene-pathway associations across multiple disease categories using GWAS data and pathway gene sets.
It identifies genes associated with each disease category, maps them to biological pathways, and finds pathways shared
between pairs of disease categories. The results are output to separate files for male and female datasets.

Functions:
----------
filter_genes_by_category(gwas_data)
    Reads a GWAS data file and filters genes based on their association with predefined disease categories.
    Returns a dictionary mapping each category to a set of associated genes.

find_pathways_for_categories(filtered_genes, pathway_data)
    Maps filtered genes for each category to pathways from a pathway gene set file.
    Returns a dictionary mapping each category to a dictionary of genes and their associated pathways,
    and a dictionary mapping pathway names to their total gene counts.

write_into_file(pathways_categorized, filename)
    Writes the dictionary including each pathway that was linked to a gene into a text file.

find_common_pathways_for_pairs(categorized_pathways, pairs)
    Identifies pathways that are common between pairs of disease categories.
    Returns a dictionary mapping each pathway to the categories and genes involved.

find_comparison_pairs(categories)
    Generates all unique pairs of disease categories for comparison.
    Returns a list of tuples, each containing a pair of categories.

create_output_files(output_file_names, categories)
    Creates or resets output files for storing results, and writes the header row.

write_output_file(common_pathways, sex, output_file_names, categories, total_genes_for_pathway)
    Writes the common pathway analysis results to the appropriate output file (male or female).
    Each row contains pathway name, category involvement, number of categories in common, total genes, and relative genes.

"""

import csv
from collections import defaultdict
import os

# Update categories here if changes occur
categories = ["cancer", "chf", "copd", "diabetes", "dementia", "mi", "stroke"]

# Input file names, CHANGE HERE
# Male data first, Female data second
pathway_data_file = 'msigdb.v2024.1.Hs.symbols.gmt'
#pathway_data_file = 'c8.all.v2024.1.Hs.symbols.gmt' #Partial data
input_file_names = ["male_gwas_ukb.txt", "female_gwas_ukb.txt"]
pathways_to_genes_file_names = ["male_pathways_to_genes_ukb.txt", "female_pathways_to_genes_ukb.txt"]
output_file_names = ['male_common_pathways_ukb.txt', 'female_common_pathways_ukb.txt']

# Other global variables
total_genes_for_pathway = set()
total_genes = -1


def filter_genes_by_category(gwas_data):
    # Filter genes based on the category
    # Dictionary with categories as keys and lists of genes as values
    category_genes = defaultdict(set)
    # Read the GWAS data file
    with open(gwas_data, newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        categories_from_file = reader.fieldnames[1:-1] # First column is 'gene' and the last column is 'incommon'
        # Iterate through each row in the file
        for row in reader:
            gene = row['gene']
            # Check each category to see if the gene is involved
            for category in categories_from_file:
                if row[category] == '1':  # check if gene is involved in the category
                    category_genes[category].add(gene)
    category_genes = dict(sorted(category_genes.items()))
    
    #returns {category : [genes]}
    return category_genes

def find_pathways_for_categories(filtered_genes, pathway_data):
    # Find pathways that include genes associated with each given category (from filtered_genes)
    # Read the pathway data file
    with open(pathway_data, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        # Read the pathway data into a list
        pathway_data = list(reader)
    # Find pathways that include genes associated with each given category
    # Initialize dictionary with the categories as keys
    pathways_categorized = {category: {} for category in categories}
    # Check if any gene in the pathway is in the filtered genes
    for category, genes in filtered_genes.items():
        for row in pathway_data:
            pathway_name = row[0] # Pathway that includes any genes that are found in GWAS data (if value = 1 in category)
            pathway_genes = set(row[2:]) # All genes in that pathway
            total_genes_for_pathway.add((pathway_name, len(pathway_genes)))
            if genes & pathway_genes:  # intersection = gene that has been associated with a pathway that belongs to a category. “Which gene names appear in both sets, exactly the same, character by character?"
                for gene in genes & pathway_genes:
                    if gene not in pathways_categorized[category]:
                        pathways_categorized[category][gene] = []
                    pathways_categorized[category][gene].append(pathway_name)
    pathways_categorized = dict(sorted(pathways_categorized.items()))

    # Returns {category: {gene:[pathways]}}
    return pathways_categorized

def write_into_file(pathways_categorized, filename):
    with open(filename, 'w') as f:
        for category, genes in pathways_categorized.items():
            f.write(f"Category: {category}\n")
            for gene, pathways in genes.items():
                f.write(f"  Gene: {gene}\n")
                for pathway in pathways:
                    f.write(f"    - Pathway: {pathway}\n")
            f.write("\n")

def find_common_pathways_for_pairs(categorized_pathways, pairs):
    global_common_pathways = defaultdict(lambda: defaultdict(list))
    # Iterate through all pairs
    for pair in pairs:
        category1 = pair[0]
        category2 = pair[1]

        pathways_with_genes_category1 = defaultdict(set)
        pathways_with_genes_category2 = defaultdict(set)

        for gene, pathways in categorized_pathways[category1].items():
            for pathway in pathways:
                pathways_with_genes_category1[pathway].add(gene)

        for gene, pathways in categorized_pathways[category2].items():
            for pathway in pathways:
                pathways_with_genes_category2[pathway].add(gene)

        # Find common pathways between genes by categories
        # Iterate through the common genes
        common_pathways = set(pathways_with_genes_category1.keys()) & set(pathways_with_genes_category2.keys())

        for pathway in sorted(common_pathways):
            global_common_pathways[pathway][category1].extend(pathways_with_genes_category1[pathway])
            global_common_pathways[pathway][category2].extend(pathways_with_genes_category2[pathway])
    # Returns {pathway_name : {category : [genes]}}
    return global_common_pathways

def find_comparison_pairs():
    pairs = []
    index1 = 0
    index2 = 1
    i = len(categories)-1
    for category in range(len(categories)):
        category1 = categories[index1]
        for j in range(i):
            category2 = categories[index2]
            pair = (category1, category2)
            pairs.append(pair)
            index2 += 1
        i -= 1
        index1 += 1
        index2 = index1 + 1

    return pairs

def create_output_files():
    for file in output_file_names:
        if os.path.isfile(file):
            os.remove(file)
            f = open(file, "x")
        # First row
        with open(file, "a") as f:
            f.write("pathway")
            f.write("\t")
            for i in range(len(categories)):
                f.write(f"{str(categories[i])}\t")
            f.write("incommon\ttotal_genes\trelative_genes")
            f.write("\n")

def write_output_file(common_pathways, sex):
    if sex == 'male':
        file_name = output_file_names[0]
    else:
        file_name = output_file_names[1]

    # Common pathways dictionary {pathway_name : {category : [genes]}}
    with open(file_name, "a") as f:
        # Iterate over each pathway that was found in common between category pairs
        for pathway, category_dict in common_pathways.items():
            row = [pathway] # First column is pathway name
            incommon = 0 # Total number of categories a pathway is involved in
            # Check presence of the category in each pathway
            for category in categories:
                if category in category_dict:
                    row.append("1")
                    incommon += 1
                    relative_genes = sum(len(genes) for genes in category_dict.values()) # Number of genes that were used to map pathway data
                else:
                    row.append("0")
            for name, count in total_genes_for_pathway:
                if name == pathway:
                    total_genes = count

            row.append(str(incommon))
            row.append(str(total_genes))
            row.append(str(relative_genes))
            f.write("\t".join(row) + "\n")


# MAIN SCRIPT

create_output_files()

# Genes associated with each category
# Returns {category : [genes]}
category_genes_male = filter_genes_by_category(input_file_names[0])
category_genes_female = filter_genes_by_category(input_file_names[1])

# Pathways associated with each category
# Returns {category: {gene:[pathways]}}
pathways_categorized_male = find_pathways_for_categories(category_genes_male, pathway_data_file)
pathways_categorized_female = find_pathways_for_categories(category_genes_female, pathway_data_file)
write_into_file(pathways_categorized_male, pathways_to_genes_file_names[0])
write_into_file(pathways_categorized_female, pathways_to_genes_file_names[1])

# Find all possible pairs from categories
# Returns [pairs]
category_pairs_male = find_comparison_pairs()
category_pairs_female = find_comparison_pairs()

# Find common pathways
# Returns {pathway_name : {category : [genes]}}
common_pathways_male = find_common_pathways_for_pairs(pathways_categorized_male, category_pairs_male)
common_pathways_female = find_common_pathways_for_pairs(pathways_categorized_female, category_pairs_female)

# Write pathways to separate files
write_output_file(common_pathways_male, 'male')
write_output_file(common_pathways_female, 'female')
