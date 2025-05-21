# ------------------------------------------------------------------------------
# Script: combine_gwas_summaries.py
# Purpose: Combine gProfiler files including rsids and gene names into one file
#          to be used in further analysis. Maps the rsids to genes and categories
#          in case needed later on.
#
# Instructions:
#   - The input file should contain a column named "id" and a column named "gene_names"
#   - Change names of the files accordingly in "female_files", "male_files", 
#     "output_filenames" and "rsid_output_files"
#
# Output:
#   - Tab separated file consisting of combined GWAS data files with columns:
#     gene, cancer, chf, copd, diabetes, dementia, mi, stroke, incommon
#   - RSID to gene mapping file
# ------------------------------------------------------------------------------

import pandas as pd

female_files = {
    "cancer": "female_cancer_gprofiler.csv",
    "chf": "female_chf_gprofiler.csv",
    "copd": "female_copd_gprofiler.csv",
    "diabetes": "female_diabetes_gprofiler.csv",
    "dementia": "female_dementia_gprofiler.csv",
    "mi": "female_mi_gprofiler.csv",
    "stroke": "female_stroke_gprofiler.csv"
}

male_files = {
    "cancer": "male_cancer_gprofiler.csv",
    "chf": "male_chf_gprofiler.csv",
    "copd": "male_copd_gprofiler.csv",
    "diabetes": "male_diabetes_gprofiler.csv",
    "dementia": "male_dementia_gprofiler.csv",
    "mi": "male_mi_gprofiler.csv",
    "stroke": "male_stroke_gprofiler.csv"
}

#NOTE order
output_file_names = ["female_gwas_ukb.txt","male_gwas_ukb.txt"]
rsid_output_files = ["rsid_to_gene_mapping_female.csv", "rsid_to_gene_mapping_male.csv"]
file_sets = [female_files, male_files]

def rsid_mapping(file_dict, output_file):
    rsid_gene_mapping = []

    for trait, file_path in file_dict.items():
        df = pd.read_csv(file_path)
        for _, row in df.iterrows():
            rsid = row['id']
            if pd.isna(row['gene_names']):
                continue
            genes = [g.strip() for g in row['gene_names'].split(",") if g.strip()]
            rsid_gene_mapping.extend((rsid, gene, trait) for gene in genes)

    df_out = pd.DataFrame(rsid_gene_mapping, columns=["rsid", "gene", "trait"]).drop_duplicates()
    df_out.to_csv(output_file, index=False)

# Function to extract gene names from the CVS files
def extract_genes(file):
    data = pd.read_csv(file)
    
    gene_list = []

    for entry in data['gene_names']:
        if pd.isna(entry):
            continue  # skip if the value is NaN
        if entry.strip() == "":
            continue  # skip if the value is an empty string

        # Split by comma and loop through each gene
        split_genes = entry.split(",")
        for gene in split_genes:
            clean_gene = gene.strip()
            if clean_gene not in gene_list:  # Duplicates away
                gene_list.append(clean_gene)

    return sorted(gene_list)

# Keep track of rsids
for i in range(len(rsid_output_files)):
    files = file_sets[i]
    output_file = output_file_names[i]
    rsid_mapping(file_sets[i], rsid_output_files[i])

# Apply to each file
gene_lists = {}

# Produce file
for i in range(len(output_file_names)):
    files = file_sets[i]
    output_file = output_file_names[i]

    # Extract genes by category
    gene_lists = {}
    for category, path in files.items():
        gene_lists[category] = extract_genes(path)

    # Get unique genes
    all_genes = []
    for gene_list in gene_lists.values():
        for gene in gene_list:
            if gene not in all_genes:
                all_genes.append(gene)
    all_genes.sort()

    # Write output file
    with open(output_file, "w") as f:
        header = ["gene"] + list(files.keys())
        f.write("\t".join(header) + "\t" + "incommon"+ "\n")
        for gene in all_genes:
            row = [gene]
            incommon = 0
            for category in files.keys():
                if gene in gene_lists[category]:
                    row.append("1")
                    incommon += 1
                else:
                    row.append("0")
            f.write("\t".join(row) + "\t")
            f.write(f"{incommon}\n")

