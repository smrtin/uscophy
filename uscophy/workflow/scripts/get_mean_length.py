from Bio import SeqIO
import sys
import numpy as np
from pathlib import Path


def main(INFILES,OUT_FILE ):

    f = open(OUT_FILE, "w")

    for FileIn in INFILES:
        #FileIn = sys.argv[1]

        handle = open(FileIn, 'r')
        sequence_lengths = {}
        SeqRecords = SeqIO.parse(handle, 'fasta')
        for record in SeqRecords:   #loop through each fasta entry
            length = len(record.seq.replace('-',''))    #get sequence length
            sequence_lengths[record.id] = length

        busco_id = Path(FileIn).stem

        #access dictionary outside of loop
        length_list=list(sequence_lengths.values())
        mean = int(np.mean(length_list)*.95)
        #standard_deviation = int(np.std(length_list)*5)

        value = int(np.var(length_list))

        f.write('{}\t0\t{}\t{}\n'.format(busco_id,value,mean))

    f.close()

if __name__ == "__main__":

    INFILES=list(snakemake.input )
    OUT_FILE=str(snakemake.output) #output needs to be named in rule!

    print(type(INFILES))
    print(type(OUT_FILE))

    main(INFILES,OUT_FILE )
