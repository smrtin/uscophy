from Bio import SeqIO
from pathlib import Path

def main(input_file, output_file):
    """
    Renames FASTA headers in an alignment file.
    New headers are formatted as: <input_file_basename>_<counter>
    Args:
        input_file (str): Path to the input FASTA file.
        output_file (str): Path to the output FASTA file with renamed headers.
    """
    # Open input and output files 
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Parse FASTA records from the input file
        seq_records = SeqIO.parse(infile, 'fasta')

        # Use the base name (without extension) of the input file as the header prefix
        file_prefix = Path(input_file).stem

        # Initialize a counter for unique header suffixes
        counter = 0

        # Loop through each FASTA record and write with new header
        for record in seq_records:
            # Create new header: e.g., >busco_id_0
            new_header = f'>{file_prefix}_{counter}'
            print(new_header, file=outfile)
            print(str(record.seq), file=outfile)
            counter += 1

if __name__ == "__main__":
    # Example usage with Snakemake input/output variables
    input_file = str(snakemake.input)
    output_file = str(snakemake.output)  # Output file name should be set in the Snakemake rule
    main(input_file, output_file)
