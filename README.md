# Northwestern Atlantic ocean microbial community difference in small and large size fractions characterised through amplicon sequence

We determined the difference in community diversity between two size fractions of water samples from a single site, and the possbility of merging two size fractions to a comparable unit to samples that were not size fractionated.

### Structure of the starting repository

    
### Datasets
The data is publically available in NCBI under project accession number [PRJNA785606](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA785606) and the ISA-tab formatted metadata is available as SraTableRun.csv.

### Generate the data
Run ```generate_data.sh``` from the cloned repo from the bash_scripts folder. This file include all the code used to process the raw FASTQ files. Then, open the ```bb2_sf_microbial.ipynb``` jupyter notebook to generate the figures and statistics which depend on ```bb2022_functions.py```. The R code was used to generate the PCoA plots, and the PERMANOVA analyses.

### Final repository
Resulting files are saved in the outputs folder.
