#!/usr/bin/env python3

import click
import pandas as pd
import os
import shutil
from zipfile import ZipFile, ZIP_DEFLATED


@click.command()
@click.argument('zipfilenames', nargs=-1, type=click.Path(exists=True))
@click.option('--output-dir', '-o', default="Reference_taxon_selection", show_default=True,
              help="Output directory for extracted files")
@click.option('--tax-file', '-t', default=None, metavar='FILE',
              help="Optional filename to write taxon ID list; if not set, taxon IDs are not written")
@click.option('--extract-gff', is_flag=True, default=False, show_default=True,
              help="Also extract GFF3 files")
@click.option('--generate-category', is_flag=True, default=False, show_default=True,
              help="Generate two-column Sample-Category metadata CSV file")

              
def main(zipfilenames, output_dir, tax_file, extract_gff, generate_category):
    """
    Extract genome FASTA files (and optionally GFF3) from NCBI dataset zip archives.
    Optionally generates a taxon ID list file (if --tax-file is specified)
    and a Sample-Category metadata CSV file.
    Creates a zipped archive of the output directory.
    """

    if not zipfilenames:
        click.echo("Error: Please provide at least one input zip filename.")
        return

    base_path = "ncbi_dataset/data"
    tmp_dir = './tmp'

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)

    tax_out = None
    if tax_file:
        tax_out = open(tax_file, "w", buffering=1)

    sample_meta = []

    for zip_i, zipfilename in enumerate(zipfilenames, 1):
        click.echo(f"Processing ZIP file {zip_i}/{len(zipfilenames)}: {zipfilename}")
        with ZipFile(zipfilename, 'r') as zip_ref:
            jsonl_path_in_zip = os.path.join(base_path, 'assembly_data_report.jsonl')
            jsonl_tmp_path = os.path.join(tmp_dir, base_path, 'assembly_data_report.jsonl')

            os.makedirs(os.path.dirname(jsonl_tmp_path), exist_ok=True)

            try:
                zip_ref.extract(jsonl_path_in_zip, tmp_dir)
            except KeyError:
                click.echo(f"Warning: {jsonl_path_in_zip} not found in {zipfilename}, skipping this archive.")
                continue

            data = pd.read_json(jsonl_tmp_path, lines=True)

            for sample_i, row in enumerate(data.itertuples(), 1):
                accession = row.accession
                assemblyName = row.assemblyInfo['assemblyName']
                organismName_orig = row.organism['organismName']
                organismTaxID = row.organism['taxId']

                sample_name = organismName_orig.replace(' ', '_')

                if tax_out:
                    tax_out.write(f"{organismTaxID}\n")

                if generate_category:
                    sample_meta.append({"Sample": sample_name, "Category": "NCBI-Genome"})

                assembly_file = f"{accession}_{assemblyName.replace(' ', '_').replace('+', '_')}_genomic.fna"
                genome_path_in_zip = os.path.join(base_path, accession, assembly_file)
                genome_tmp_path = os.path.join(tmp_dir, genome_path_in_zip)

                try:
                    zip_ref.extract(genome_path_in_zip, tmp_dir)
                except KeyError:
                    click.echo(f"Warning: Genome file {genome_path_in_zip} missing in {zipfilename}. Skipped sample {sample_name}.")
                    continue

                genome_out_file = f"{sample_name}.fas"
                genome_out_path = os.path.join(output_dir, genome_out_file)
                shutil.copy(genome_tmp_path, genome_out_path)
                os.remove(genome_tmp_path)

                click.echo(f"  [{sample_i}/{len(data)}] Extracted genome for sample '{sample_name}' to '{genome_out_file}'")

                if extract_gff:
                    gff_path_in_zip = os.path.join(base_path, accession, 'genomic.gff')
                    gff_tmp_path = os.path.join(tmp_dir, gff_path_in_zip)
                    try:
                        zip_ref.extract(gff_path_in_zip, tmp_dir)
                        gff_out_file = f"{sample_name}.gff3"
                        gff_out_path = os.path.join(output_dir, gff_out_file)
                        shutil.copy(gff_tmp_path, gff_out_path)
                        os.remove(gff_tmp_path)
                        click.echo(f"    Extracted GFF3 for sample '{sample_name}'")
                    except KeyError:
                        click.echo(f"    Warning: GFF3 file {gff_path_in_zip} missing for sample '{sample_name}'")

            if os.path.exists(jsonl_tmp_path):
                os.remove(jsonl_tmp_path)

    if tax_out:
        tax_out.close()

    try:
        shutil.rmtree(tmp_dir)
    except Exception:
        pass

    if generate_category and sample_meta:
        meta_df = pd.DataFrame(sample_meta)
        meta_out_file = os.path.join(output_dir, "samples_metadata.csv")
        meta_df.to_csv(meta_out_file, index=False)
        click.echo(f"Sample category metadata CSV file generated at: {meta_out_file}")


    click.echo("Extraction complete.")
    click.echo(f"Files saved in: {output_dir}")
    if tax_file:
        click.echo(f"Taxon IDs saved to: {tax_file}")


if __name__ == "__main__":
    main()
