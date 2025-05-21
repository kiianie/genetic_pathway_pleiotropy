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



### [`phi_coefficient_script.R`](phi_coefficient_script.R)

**Inputs**

- Pathway-category association file

**Outputs**

- Pairwise phi coefficient table

- Correlation plots (PNG Files) for visualization


## How to Run

**1. Prepare Input Files**
- Ensure GWAS gene-category files and pathway gene set files are present in the project directory.
- The GWAS data should include the following columns separated by tabulator: gene, cancer, chf, copd, diabetes, dementia, mi, stroke, incommon (if changes are done on the categories the pathway_pleiotropy.py and phi_coefficient.R scripts should be edited accordingly)

**2. Run pathway_pleoitropy.py**
- Edit the input and output file names to match the needs

**3. Run phi_coefficient.R**
- Input file is the common_pathway.txt file derived from Step 2.
- Run code separately for each file by changing the file name in the code
