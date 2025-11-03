from Bio import SeqIO
from pathlib import Path
import sys 

def main(INFILES,OUT_FILE ):

    fout = open(OUT_FILE, 'w')

    for FileIn in INFILES:

        handle = open(FileIn, 'r')
        header_base = Path(FileIn).stem.split('.', 1)[0]
        counter=0
        Lines = handle.readlines()

        for line in Lines:
            line.strip
            header , sequence = line.split()
            number = int(header.replace('Node','')) -1
            new_header = '>{}_{}'.format(header_base,number)
            fout.write('{}\n{}\n'.format(new_header, sequence))
    fout.close()

if __name__ == "__main__":

    INFILES=list(snakemake.input )
    OUT_FILE=str(snakemake.output) #output needs to be named in rule!

    main(INFILES,OUT_FILE )

