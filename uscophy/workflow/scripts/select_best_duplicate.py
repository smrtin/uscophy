#!/usr/bin/env python3

import sqlite3
import os
import logging
import subprocess
import tempfile
import shutil

def ensure_column_exists(db_path, column_name):
    with sqlite3.connect(db_path) as conn:
        cur = conn.cursor()
        cur.execute("PRAGMA table_info(busco_sequences);")
        existing_columns = [row[1] for row in cur.fetchall()]
        if column_name not in existing_columns:
            cur.execute(f"ALTER TABLE busco_sequences ADD COLUMN {column_name} INTEGER DEFAULT 0;")
            conn.commit()

def copy_database(in_path, out_path):
    shutil.copy2(in_path, out_path)

def extract_and_update_duplicates(database, diamond_database, threads, tmp_base_dir="local_tmp"):
    ensure_column_exists(database, "best_duplicate")

    os.makedirs(tmp_base_dir, exist_ok=True)

    with sqlite3.connect(database) as db_connection:
        cursor = db_connection.cursor()
        cursor.execute(
            "SELECT b.Busco_id, b.Sequence, b.Species, b.Header FROM busco_sequences b "
            "INNER JOIN full_table f ON b.Busco_id = f.Busco_id AND b.Species = f.Species "
            "WHERE Status='Duplicated' AND b.Sequence_type = 'amino' "
            "ORDER BY b.Busco_id, b.Species;"
        )
        rows = cursor.fetchall()
        grouped = {}
        for busco_id, sequence, species, header in rows:
            grouped.setdefault((busco_id, species), []).append((sequence, header))

        for (busco_id, species), seqs in grouped.items():
            if len(seqs) > 1:
                print(f"Checking for duplicated sequences in {busco_id} ({species})")
                with tempfile.TemporaryDirectory(dir=tmp_base_dir) as tmpdir:
                    fasta_path = os.path.join(tmpdir, f"{busco_id}_{species}_query.fas")
                    result_path = os.path.join(tmpdir, f"{busco_id}_{species}_dmnd.m8")
                    with open(fasta_path, "w") as fa:
                        for seq, header in seqs:
                            fa.write(f">{header}\n{seq}\n")

                    subprocess.run([
                        "diamond", "blastp",
                        "-q", fasta_path,
                        "-d", diamond_database,
                        "-o", result_path,
                        "--threads", str(threads) ,
                        "--fast", "--max-target-seqs", "1", "--outfmt", "6" 
                    ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    # maybe also use sensitive mode

                    with open(result_path) as f:
                        best_score = -1
                        best_header = None
                        for line in f:
                            parts = line.strip().split('\t')
                            qname, score = parts[0], float(parts[11])
                            
                            if score > best_score:
                                best_score = score
                                best_header = qname

                    if best_header is not None:
                        cursor.execute(
                            "UPDATE busco_sequences SET best_duplicate=1 WHERE Header=?",
                            (best_header,))
                        for _, header in seqs:
                            if header != best_header:
                                cursor.execute(
                                    "UPDATE busco_sequences SET best_duplicate=0 WHERE Header=?",
                                    (header,))
        db_connection.commit()

if __name__ == '__main__':
    orig_database = snakemake.input.busco_db
    diamond_database = snakemake.input.diamond_db

    output_db = snakemake.output[0]
    logfile = snakemake.log[0]
    threads = snakemake.threads

    logging.basicConfig(filename=logfile,
                       level=logging.INFO,
                       format="%(asctime)s %(levelname)s:%(message)s")

    copy_database(orig_database, output_db)
    extract_and_update_duplicates(output_db, diamond_database, threads, tmp_base_dir="local_tmp")
