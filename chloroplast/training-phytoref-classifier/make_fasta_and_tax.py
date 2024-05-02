#generate the two files necessary to train a classifier on phytoref sequencing data
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import SeqIO
#import fasta file
fa = SeqIO.parse('chloroplast/training-phytoref-classifier/euk_cyano.fasta' ,
                 "fasta")
#create empty list for original_ids
original_ids = []
#make a list of the original ids from PhytoRef
for record in fa: #a SeqRecord has the accession as record.id, usually.
    original_ids.append(record.id)

#switch the list to df
df = pd.DataFrame({'originalid':original_ids})
#split the name into its id and taxonomy
df[['newid','taxonomy']] = df['originalid'].str.split('|', 1, expand=True)

#make dic for new id assignment
edit_id = pd.Series(df.newid.values,df.originalid.values).to_dict()

#use dic to edit the original fasta file into a corrected one
with open('chloroplast/training-phytoref-classifier/euk_cyano.fasta') as original, open('chloroplast/training-phytoref-classifier/euk_cyano_idonly.fasta', 'w') as corrected:
    for seq_record in SeqIO.parse(original, 'fasta'):
        if seq_record.id in edit_id:
            seq_record.id = seq_record.description = edit_id[seq_record.id]
        SeqIO.write(seq_record, corrected, 'fasta')
#export the taxonomy file
df[['newid', 'taxonomy']].copy().to_csv('chloroplast/training-phytoref-classifier/euk_cyano_taxo.tsv', sep='\t', index=False, header=False)
