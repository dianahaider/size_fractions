
#import the sequnces from cyano/chloroplast from whole dataset
qiime tools import \
  --input-path selected_dna.fasta \
  --output-path rep_sequences.qza \
  --type 'FeatureData[Sequence]'

#make a folder to train the classifier
mkdir training-phytoref-classifier
cd training-phytoref-classifier

#import the fasta and taxonomy tsv for phytoref

wget http://phytoref.sb-roscoff.fr/static/downloads/Cyanobacteria_p98.fasta
wget http://phytoref.sb-roscoff.fr/static/downloads/PhytoRef_with_taxonomy.fasta

cat * >> euk_cyano.fasta

#there was one occurence of 10Xs, I replaced them with 10Ns
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path euk_cyano_idonly.fasta \
  --output-path euk_cyano.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path euk_cyano_taxo.tsv \
  --output-path ref-taxonomy.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads euk_cyano.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier classifier.qza

qiime feature-classifier classify-sklearn \
  --i-classifier training-phytoref-classifier/classifier.qza \
  --i-reads rep_sequences.qza \
  --o-classification taxonomy_chloro.qza

#extract only the taxonomy
for zipfile in taxonomy_chloro.qza; do unzip -j "$zipfile" '*.tsv' -x '*/*/*/*'; done
