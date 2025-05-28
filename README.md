# Pathway Pleiotropy Analysis

This project analyzes gene-pathway associations across multiple disease categories using GWAS data and pathway gene sets. It identifies genes associated with each disease category, maps them to biological pathways, and finds pathways shared between pairs of disease categories. The results are output to separate files for male and female datasets.

# Used software
Python version 3.11.9.

R version 4.5.0.

## Main Scripts

### [`pathway_pleiotropy.py`](pathway_pleiotropy.py)

**Inputs**

- GWAS datasets (divided by sex)

- Pathway gene set

**Outputs**

- Pathway-category association files divided by sex
- (Secondary output : pathways mapped to genes by trait)


### [`phi_coefficient_script.R`](phi_coefficient_script.R)

**Inputs**

- Pathway-category association file

**Outputs**

- Tab separated file with pairwise phi coefficients

- Phi correlation matrix and correlation plot

### [`phi_comparison_table_script.R`](phi_comparison_table_script.R)

**Inputs**

- Two pairwise phi coefficient tables

**Outputs**

- Tab separated file with phi values on from both files

### [`combine_gwas_summaries.py`](combine_gwas_summaries.py)

**Inputs**

- gProfiler files of each category

**Outputs**

-  Tab separated file consisting of combined GWAS data files with columns: gene, cancer, chf, copd, diabetes, dementia, mi, stroke, incommon
-  RSID to gene mapping file

## How to Run

**1. Prepare Input Files**
- Ensure GWAS gene-category files and pathway gene set files are present in the project directory.
- The GWAS data should include the following columns separated by tabulator: gene, cancer, chf, copd, diabetes, dementia, mi, stroke, incommon (if changes are done on the categories the pathway_pleiotropy.py and phi_coefficient.R scripts should be edited accordingly)
- **If needed, gProfiler files derived from GWAS summaries can be combined with combine_gwas_summaries.py**

**2. Run pathway_pleoitropy.py**
- Edit the input and output file names
- Run code

**3. Run phi_coefficient.R**
- Input file is the common_pathway.txt file derived from Step 2.
- Edit the input and output file names
- Run code separately for each file by changing the file name in the code

**Optional: 4. Run phi_comparison_table_script.R**
- Input files are the output files from Step 3.
- Edit input and output file names
- Run code
