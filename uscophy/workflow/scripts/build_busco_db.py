#!/usr/bin/env python3

import os
import sqlite3
import csv 
import pandas as pd
import glob
from Bio import SeqIO
import numpy as np


def build_empty_db(data_path, filename ):
    db = sqlite3.connect(os.path.join(data_path,filename))

    db.execute('CREATE TABLE "BuscoSammlung" (\
	"Species_name"	TEXT,\
	"Busco_id"	TEXT,\
	"Status"	TEXT,\
	"Sequence"	TEXT,\
	"Gene_Start"	TEXT,\
	"Gene_End"	INTEGER,\
	"Strand"	TEXT,\
	"Score"	INTEGER,\
	"Length"	INTEGER,\
	"OrthoDB_url"	TEXT,\
	"Description"	TEXT  \
    );')
    db.close()


def load_csv(input_table):
    with open(input_table,'r') as f_in:
        next(f_in)
        next(f_in)
        df = pd.read_csv(f_in, delimiter='\t')
        df.columns = [c.replace('# ', '') for c in df]
        df.columns = [c.replace(' ', '_') for c in df]

    return df

def read_in_single_copy(busco_path, sequence_type):

    sequence_type_extension = {
        'amino': '.faa',
        'nucleotide': '.fna'
    }

    sequence_list = []

    for file in glob.iglob(os.path.join(busco_path,'busco_sequences/single_copy_busco_sequences/*{}'.format(sequence_type_extension[sequence_type]))):
        base=os.path.basename(file)
        busco_id = os.path.splitext(base)[0]
        
        FastaFile = open(file, 'r')
        for seqs in SeqIO.parse(FastaFile, 'fasta'):
            name = seqs.id
            sequence=seqs.seq

            sequence_list.append((busco_id,sequence_type,name,str(sequence)))

    return sequence_list

def read_in_fragmented(busco_path, sequence_type):

    sequence_type_extension = {
        'amino': '.faa',
        'nucleotide': '.fna'
    }

    sequence_list = []

    for file in glob.iglob(os.path.join(busco_path,'busco_sequences/fragmented_busco_sequences/*{}'.format(sequence_type_extension[sequence_type]))):
        base=os.path.basename(file)
        busco_id = os.path.splitext(base)[0]
        
        FastaFile = open(file, 'r')
        for seqs in SeqIO.parse(FastaFile, 'fasta'):
            name = seqs.id
            sequence=seqs.seq

            sequence_list.append((busco_id,sequence_type,name,str(sequence)))

    return sequence_list

def read_in_duplicated(busco_path, sequence_type):

    sequence_type_extension = {
        'amino': '.faa',
        'nucleotide': '.fna'
    }

    sequence_list = []

    for file in glob.iglob(os.path.join(busco_path,'busco_sequences/multi_copy_busco_sequences/*{}'.format(sequence_type_extension[sequence_type]))):
        base=os.path.basename(file)
        busco_id = os.path.splitext(base)[0]
        
        FastaFile = open(file, 'r')
        for seqs in SeqIO.parse(FastaFile, 'fasta'):
            name = seqs.id
            sequence=seqs.seq

            sequence_list.append((busco_id,sequence_type,name,str(sequence)))

    return sequence_list



def main(input_path,output_path,sample):



    filename = 'busco_database.sqlite3'

    os.makedirs(output_path, exist_ok=True)
    databasename = os.path.join(output_path,filename)

    if os.path.exists(databasename): 
        os.remove(databasename)

    db = sqlite3.connect(databasename)

    ## load busco results table
    busco_table_path = os.path.join(input_path,'full_table.tsv')
    data = load_csv(busco_table_path)
    data['Species'] = sample
    
    ## Some columns are missing in the transcriptome output

    if 'Gene_Start' not in data.columns:
        data["Gene_Start"] = np.nan
    if 'Gene_End' not in data.columns:
        data["Gene_End"] = np.nan
    if 'Strand' not in data.columns:
        data["Strand"] = np.nan    

    busco_table_name = 'full_table'
    data.to_sql(busco_table_name, db, if_exists='replace')

    ## load busco results sequences 
    sequence_list = []

    sequence_list.extend(read_in_single_copy(input_path, 'amino' ))
    sequence_list.extend(read_in_single_copy(input_path, 'nucleotide'))
    sequence_list.extend(read_in_fragmented(input_path, 'amino' ))
    sequence_list.extend(read_in_fragmented(input_path, 'nucleotide'))
    sequence_list.extend(read_in_duplicated(input_path, 'amino' ))
    sequence_list.extend(read_in_duplicated(input_path, 'nucleotide'))

    sequence_data = pd.DataFrame(sequence_list, columns=('Busco_id', 'Sequence_type', 'Header', 'Sequence'))
    sequence_data['Species'] = sample

    busco_table_name = 'busco_sequences'
    sequence_data.to_sql(busco_table_name, db, if_exists='replace')


if __name__ == "__main__":

    IN_PATH=str(snakemake.params.input_dir)
    OUT_PATH=str(snakemake.params.input_dir) #output needs to be named in rule!
    SAMPLE=str(snakemake.wildcards.sample)

    main(IN_PATH, OUT_PATH,SAMPLE )
