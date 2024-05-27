#1.0 Inspect read quality
mkdir fastqc_out
fastqc -t $NCORES 00-raw/*.fastq.gz -o fastqc_out

cd fastqc_out
multiqc .
cd ..

#2.0
conda activate cutadaptenv

source ~/Documents/escuela/phd/eASV-pipeline-for-515Y-926R/qiime2-2022.2-DADA2-SILVA138.1-PR2_4.14.0/00-trimming-sorting-scripts/00-run-cutadapt.sh















#2.0 Download qiime environment and activate
#follow download instructions on https://docs.qiime2.org/2023.5/install/native/
source activate qiime2-2023.5 

mkdir reads_qza
    
qiime tools import \
   --type SampleData[PairedEndSequencesWithQuality] \
   --input-path raw_data/ \
   --output-path reads_qza/reads.qza \
   --input-format CasavaOneEightSingleLanePerSampleDirFmt

qiime cutadapt trim-paired \
   --i-demultiplexed-sequences reads_qza/reads.qza \
   --p-cores 1 \
   --p-front-f GTGYCAGCMGCCGCGGTAA \
   --p-front-r CCGYCAATTYMTTTRAGTTT \
   --p-discard-untrimmed \
   --p-no-indels \
   --o-trimmed-sequences reads_qza/reads_trimmed.qza


#from local: 
ssh -Y dhaider@narval.computecanada.ca
(to ssh my pssword is 9043..)
scp -r 2017 dhaider@narval.computecanada.ca:~/projects/def-rbeiko/dhaider/bb_data/split_pipeline

#**move to compute canada**
#from remote:
singularity exec -B $PWD:/home/projects/def-rbeiko/dhaider/bb_data -B 2017:/outputs \
  -B 2017:/inputs qiime2-2023.2.sif \
  qiime dada2 denoise-paired --i-demultiplexed-seqs /inputs/reads_trimmed.qza --p-trunc-len-f 280 --p-trunc-len-r 280 --output-dir /outputs/dada2_output --verbose

wget https://data.qiime2.org/2023.5/common/silva-138-99-nb-classifier.qza

singularity --mem=16GB exec -B $PWD:/home/projects/def-rbeiko/dhaider/bb_data -B 2014:/outputs \
  -B 2014/dada2_output:/inputs -B 2014:/class qiime2-2023.2.sif \
  qiime feature-classifier classify-sklearn --i-reads /inputs/representative_sequences.qza --i-classifier /class/silva-138-99-nb-classifier.qza --output-dir /outputs/taxa --verbose


#from local
qiime feature-table summarize --i-table dada2_output_270210/table.qza --o-visualization dada2_output_270210/dd2270210_table_summary.qzv

qiime feature-table filter-features --i-table dada2_output_270210/table.qza --p-min-frequency 34 --p-min-samples 1 --o-filtered-table dada2_output_270210/table_filt.qza

qiime taxa filter-table --i-table dada2_output_270210/table_filt.qza --i-taxonomy taxa_270210/classification.qza --p-exclude mitochondria,chloroplast --o-filtered-table dada2_output_270210/table_filt_contam.qza

qiime feature-table summarize --i-table dada2_output_270210/table_filt_contam.qza --o-visualization dada2_output_270210/table_filt_contam_summary.qzv


qiime diversity alpha-rarefaction --i-table dada2_output_270210/table_filt_contam.qza --p-max-depth 34 --p-steps 20 --p-metrics 'observed_features' --o-visualization rarefaction_curves_test_270210.qzv

qiime feature-table filter-samples --i-table dada2_output_270210/table_filt_contam.qza --p-min-frequency 7500 --o-filtered-table dada2_output_270210/table_filt_min.qza


qiime feature-table filter-seqs --i-data dada2_output_270210/representative_sequences.qza --i-table dada2_output_270210/table_filt_contam.qza --o-filtered-data dada2_output_270210/rep_seqs_filt_contam_final.qza


qiime fragment-insertion sepp --i-representative-sequences dada2_output_270210/rep_seqs_filt_contam_final.qza --i-reference-database sepp-refs-gg-13-8.qza --o-tree asvs-tree.qza --o-placements insertion-placements.qza


qiime diversity core-metrics-phylogenetic --i-table dada2_output_270210/table_filt_contam.qza --i-phylogeny asvs-tree.qza --p-sampling-depth 8000 --m-metadata-file METADATA_2.tsv --p-n-jobs-or-threads 1 --output-dir diversity --verbose


qiime composition ancom --i-table dada2_output_270210/table_filt_contam_pseudocount.qza --m-metadata-file METADATA_2.tsv --m-metadata-column 'Depth code' --output-dir ancom_output

#export biom table
qiime tools export    --input-path dada2_output/table_filt_contam.qza    --output-path dada2_output_exported

biom convert -i feature-table.biom -o feature-table.tsv --to-tsv
