#!/usr/bin/env python3

import os, sys
import subprocess
import click
import re
from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn
from snakemake.io import glob_wildcards
import shlex


def get_snakefile(file="workflow/Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    # print(sf)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf

def get_snakefile_build(file="workflow/Snakefile_build"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    # print(sf)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file for building busco sets; tried %s" % sf)
    return sf

def get_snakefile_assemble(file="workflow/Snakefile_assemble"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    # print(sf)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file for target assembly; tried %s" % sf)
    return sf

@click.group()
def cli():
    pass

## build busco set

@cli.command(
    "build"
)
@click.option(
    "-t", 
    "--threads",  
    type=int,
    help="provide number of threads / cores / jobs ", 
    default=5
)
@click.option(
    "-d", 
    "--dry",  
    is_flag=True,
    default=False,
    show_default=True,
    help="Test execution.", 
)
@click.option(
    "-i", 
    "--input",    
    help="directory that contains all alignment files with a .fas ending [ default: input ]",
    default='input'
)
@click.option(
    "-o", 
    "--output",    
    help="directory that contains the busco sets [ default: ./ ]",
    default='lineages'
)
@click.option(
    "-n",
    "--set_name",
    help="name your busco set [default: custom_set_odb10 ]",
    default='custom_set_odb10'
)
@click.option(
    "-s",
    "--snakemake",  
    help="additional snakemake options and command", 
    multiple=True
)
def run_build_busco_set(threads, dry, input , output ,set_name,snakemake):

    if not set_name.endswith("_odb10"):
        set_name = set_name + '_odb10'
        print('Your setname does not end with \'_odb10\' and this might cause problems with the busco call. \n\nWe adjust it to: {} '.format(set_name))


    cmd = (
        "snakemake --nolock {dryrun} --snakefile {snakefile} --config input_dir='{input}' output_dir='{output}' set_name='{set_name}' --use-conda {snakemake}"
    ).format(
        snakefile=get_snakefile_build(),
        cores="--cores {}".format(threads) if threads is not None else "",
        dryrun="--dryrun" if dry else "",
        input=input,
        output=output, 
        set_name=set_name,
        snakemake=' '.join(snakemake) if snakemake else ""
        )  

    #print(cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:

        exit(1)




## target assembly 
@cli.command(
    "assemble"
)
@click.option(
    "-l", 
    "--lineage",  
    help="provide a busco lineage [ default: metazoa_odb10 ]", 
    default='metazoa_odb10'
)
@click.option(
    "-o", 
    "--output",    
    help="directory that contains the busco sets [ default: output ]",
    default='output'
)
@click.option(
    "--samples",    
    help="sample file with information about sample and read files",
    required=True,
)
@click.option(
    "-t", 
    "--threads",  
    type=int,
    help="provide number of threads / cores / jobs [ default: 5 ]", 
    default=5
)
@click.option(
    "-a",
    "--assembly_threads",
    type=int,
    help="provide number of threads for each assembly step [ default: 10 ]",
    default=10
)
@click.option(
    "-r",
    "--ref_genome",
    help="provide a reference genome to scaffold fragmented assembly parts", 
)
@click.option(
    "-d", 
    "--denovo",  
    is_flag=True,
    default=False,
    show_default=True,
    help="Run denovo assembly from trimmed reads. Without prefiltering reads with diamond", 
)
@click.option( 
    "--dry",  
    is_flag=True,
    default=False,
    show_default=True,
    help="Test execution.", 
)

@click.option(
    "-s",
    "--snakemake",  
    help="additional snakemake options and command", 
    multiple=True
)
def run_target_assembly(lineage, output, samples, threads,ref_genome, dry , snakemake ,denovo,assembly_threads):
    
    if os.path.exists(os.path.dirname(lineage)):
        print('is a local busco set')
        # now split into path and lineage name
        lineage_name = os.path.basename(os.path.normpath(lineage))
        
        lineage_path = os.path.dirname(os.path.normpath(lineage))
    else:
        lineage_name = lineage 
        lineage_path = os.path.join(output, 'busco/lineages')
    
    print(lineage_name) 
    print(lineage_path)

    cmd = (
        "snakemake --nolock --snakefile {snakefile} {cores} {dryrun} --config lineage='{lineage_name}' lineage_path='{lineage_path}' output_dir='{output}' {denovo} {assembly_threads} assembly_sample_info='{samples}' {ref_genome} {snakemake}"
    ).format(
        snakefile=get_snakefile_assemble(),
        lineage_name=lineage_name,
        lineage_path=lineage_path,
        output=output, 
        samples=samples,
        dryrun="--dryrun" if dry else "",
        cores="--cores {}".format(threads) if threads is not None else "",
        ref_genome="reference_genome={}".format(ref_genome) if ref_genome else "",
        denovo="denovo=True" if denovo else "",
        snakemake=' '.join(snakemake) if snakemake else "",
        assembly_threads='assembly_threads={}'.format(assembly_threads),

        )  

    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:

        exit(1)


## run workflow 

@cli.command(
    "run"
)
@click.option(
    "-g", 
    "--genomic",    
    help="directory that contains all genome files with a .fas ending [ default: input ]", 
    default='input'
)
@click.option(
    "-t", 
    "--transcriptomic",    
    help="directory that contains all transcriptomic files with a .fas ending [ default: transcriptomes ]", 
    default='transcriptomes'
)
@click.option(
    "-l", 
    "--lineage",  
    help="provide a busco lineage [ default: metazoa_odb10 ]", 
    default='metazoa_odb10'
)
@click.option(
    "-o", 
    "--output",    
    help="directory that contains the busco sets [ default: output ]",
    default='output'
)
@click.option(
    "-n", 
    "--num_taxa",  
    type=int,
    help="minimum number of taxa in alignment of the busco gene to be included in further analysis [ default: 5 ]", 
    default=5
)
@click.option(
    "-f", 
    "--frag",     
    help="also include fragmented busco genes", 
    is_flag=True,
    default=False,
    show_default=True,
)
@click.option(
    "-b", 
    "--best_duplicated",     
    help="also include best duplicated sequencs. Most similar to the ancestral variants of the busco lineage.", 
    is_flag=True,
    default=False,
    show_default=True,
)

@click.option(
    "-m", 
    "--min_genes",  
    type=float,
    help="number of minimum busco genes per sample [ default: 100 ]. When a number between 0 and 1 is provided, it will be interpreted as percentage of the total number of genes available in the busco set.", 
    default=100
)

@click.option(
    "--alignment_software", 
    help="choose an alignment software (mafft or hmmalign)",
    default='mafft',
    show_default=True, 
    type=click.Choice(['mafft', 'hmmalign'])
     )

@click.option(
    "--genetree",     
    help="also generate genetrees and reconstuct phylogeny with astral", 
    is_flag=True,
    default=False,
    show_default=True,
)
@click.option(
    '--modeltesting', 
    help='Modeltesting option for the tree reconstruction step (JTT, TEST, MFP)', 
    default='JTT', 
    show_default=True, 
    type=click.Choice(['JTT', 'TEST', 'MFP'])
)

@click.option(
    "--threads",  
    type=int,
    help="provide number of threads / cores / jobs ", 
    default=5
)
@click.option(
    "-d", 
    "--dry",  
    is_flag=True,
    default=False,
    show_default=True,
    help="Test execution.", 
)

@click.option(
    "-s",
    "--snakemake",  
    help="additional snakemake options and command", 
    multiple=True
)

@click.option(
    '--category-csv', 
    type=click.Path(exists=True), 
    default=None, 
    help='Optional: Path to sample_category.csv.'
    )

@click.option(
    '--min-taxa-per-category', 
    default=None, 
    help="Optional: Minimum taxa per category string, e.g. 'cat1:3,cat2:5'"
    )

@click.option(
    '--outgroup', 
    default=None, 
    help="Optional: Provide a Sample-ID that is used as an Outgroup in the final tree"
    )


def run_workflow(
        threads, 
        dry , 
        lineage, 
        output, 
        num_taxa,
        min_genes,
        frag , 
        best_duplicated,
        snakemake, 
        genomic,
        transcriptomic, 
        genetree, 
        alignment_software, 
        modeltesting,
        category_csv, 
        min_taxa_per_category,
        outgroup
        ):

    if os.path.exists(os.path.dirname(lineage)):
        print('is a local busco set')
        # now split into path and lineage name
        lineage_name = os.path.basename(os.path.normpath(lineage))
        
        lineage_path = os.path.dirname(os.path.normpath(lineage))
    else:
        lineage_name = lineage 
        lineage_path = os.path.join(output, 'busco/lineages')
    

    cmd = (
        "snakemake {snakemake} --latency-wait 60 --nolock {dryrun} --snakefile {snakefile} {cores} --config genomic_dir='{genomic}' genetree='{genetree}' {transcriptomic_in} lineage='{lineage_name}' lineage_path='{lineage_path}' output_dir='{output}' min_spec='{min_species}' min_genes='{min_genes}' alignment_software={alignment_software} modeltesting={modeltesting} {category_csv} {min_taxa_per_category} {outgroup} {fragments} {best_duplicated}"
    ).format(
        snakefile=get_snakefile(),
        cores="--cores {}".format(threads) if threads is not None else "",
        dryrun="--dryrun" if dry else "",
        lineage_name=lineage_name,
        lineage_path=lineage_path,
        output=output, 
        genomic=genomic,
        transcriptomic_in="transcriptomic_dir={}".format(transcriptomic)  if transcriptomic else "",
        fragments="fragmented='True'" if frag else "fragmented='False'",
        best_duplicated="duplicated='True'" if best_duplicated else "duplicated='False'" ,
        min_species=num_taxa,
        min_genes=min_genes,
        alignment_software=alignment_software,
        genetree=genetree,
        modeltesting=modeltesting,
        category_csv="sample_category_csv={}".format(category_csv) if category_csv else "",
        min_taxa_per_category="min_taxa_per_category='{}'".format(min_taxa_per_category) if min_taxa_per_category else "",
        outgroup="outgroup='{}'".format(outgroup) if outgroup else "",
        snakemake=' '.join(snakemake) if snakemake else ""
        )  
        ## --use-conda can be added via snakemake option

    
    print(shlex.split(cmd))
    print(cmd)
    #set printing color


    num_input_files=len(glob_wildcards(os.path.join(genomic , "{sample}.fas")).sample)
    num_trans_input_files=len(glob_wildcards(os.path.join(transcriptomic , "{sample}.fas")).sample)

    logo = [
        r'''  _    _                     _           ''',
        r''' | |  | |                   | |          ''',
        r''' | |  | |___  ___ ___  _ __ | |__  _   _ ''',
        r''' | |  | / __|/ __/ _ \| '_ \| '_ \| | | |''',
        r''' | |__| \__ \ (_| (_) | |_) | | | | |_| |''',
        r'''  \____/|___/\___\___/| .__/|_| |_|\__, |''',
        r'''                      | |           __/ |''',
        r'''                      |_|          |___/ '''
        ]                                                    

    


    print('\n'.join(logo))
    print( "\n - running analysis with lineage:\n\n\t{}\n".format(lineage_name) )
   
    if num_input_files == 0 and num_trans_input_files ==0 :
            print(" no files detected in input folder\n" )
            sys.exit(1)
    else:
        print(" - processing {} ".format(num_input_files) +"genomic samples\n")
        print(" - processing {} ".format(num_trans_input_files)+"transciptomic samples\n")

    try:
        with Progress(
                TextColumn("[progress.description]{task.description}"),
                BarColumn(bar_width=100),
                TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
                TimeElapsedColumn(),
                transient=True
            ) as progressbar:

                task = progressbar.add_task("Processing...", total=100)
                # Start a subprocess
                process = subprocess.Popen(
                    #cmd.split(), stderr=subprocess.PIPE, stdout=subprocess.PIPE, text=True
                    shlex.split(cmd) , stderr=subprocess.PIPE, stdout=subprocess.PIPE, text=True
                ) 
                err = 0
                while True:
                    process_out = process.stderr.readline()
                    if process_out == '' and process.poll() is not None:
                        break
                    if process_out:
                        if err > 0:
                            if "Shutting down, this" in process_out:
                                sys.exit(1)

                            print( f"{process_out.strip()}" )

                        if "Error in rule" in process_out:

                            process.kill()

                            print(f"{process_out.strip()}" )

                            err += 1
                        
                        match = re.search(r"\(\d+%\)", process_out)

                        if match:

                            percent = int(re.sub(r'\D', '', match.group()))

                            progressbar.update(task, completed=percent)

        if process.returncode < 1:

           print("\nfinished successfully\nResults were written to:\t{}\n".format(output) )

    except KeyboardInterrupt:

        # Handle the keyboard interrupt

        print("\nTerminating uscophy ...\n")

        process.kill()

        sys.exit(1)


def main():

    cli()

if __name__ == "__main__":

    main()
