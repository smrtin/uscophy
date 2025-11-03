#!/usr/bin/env python3

import sqlite3
import os
import sys
import logging

def get_column_names(cursor, table, db_prefix=None):
    """
    Fetch column names for a given table, optionally from an attached database.
    """
    if db_prefix:
        pragma_stmt = f'PRAGMA {db_prefix}.table_info("{table}")'
    else:
        pragma_stmt = f'PRAGMA table_info("{table}")'
    cursor.execute(pragma_stmt)
    return [row[1] for row in cursor.fetchall()]

def create_table_from_source(cursor, table, source_db_prefix):
    """
    Create a table in the destination database using the schema from the source database.
    """
    cursor.execute(
        f"SELECT sql FROM {source_db_prefix}.sqlite_master WHERE type='table' AND name=?",
        (table,)
    )
    create_stmt = cursor.fetchone()
    if create_stmt and create_stmt[0]:
        cursor.execute(create_stmt[0])

def merge_table(cursor, table, source_db_prefix):
    """
    Merge data from a table in the attached source database into the destination database.
    Only common columns are merged, preserving the destination column order.
    """
    dest_cols = get_column_names(cursor, table)
    src_cols = get_column_names(cursor, table, source_db_prefix)

    # Find common columns in the order of the destination table
    common_cols = [col for col in dest_cols if col in src_cols]
    if not common_cols:
        logging.info(f"Skipping table '{table}': no common columns.")
        #print(f"Skipping table '{table}': no common columns.")
        return

    cols_str = ", ".join(f'"{col}"' for col in common_cols)
    sql = (
        f'INSERT OR IGNORE INTO "{table}" ({cols_str}) '
        f'SELECT {cols_str} FROM {source_db_prefix}."{table}"'
    )
    #print(f"Merging table '{table}' with columns: {common_cols}")
    cursor.execute(sql)

def merge_databases(dest_db_path, source_db_path, tables_to_merge):
    """
    Merge specified tables from the source database into the destination database.
    """
    conn = sqlite3.connect(dest_db_path)
    cursor = conn.cursor()
    conn.execute(f"ATTACH DATABASE ? AS srcdb", (source_db_path,))
    conn.execute("BEGIN")

    for table in tables_to_merge:
        # Check if table exists in source
        cursor.execute(
            "SELECT name FROM srcdb.sqlite_master WHERE type='table' AND name=?",
            (table,)
        )
        if not cursor.fetchone():
            loggin.info(f"Table '{table}' does not exist in {source_db_path}, skipping.")
            #print(f"Table '{table}' does not exist in {source_db_path}, skipping.")
            continue

        # Create table in destination if it doesn't exist
        cursor.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
            (table,)
        )
        if not cursor.fetchone():
            logging.info(f"Creating table '{table}' in destination database.")
            #print(f"Creating table '{table}' in destination database.")
            create_table_from_source(cursor, table, "srcdb")
            conn.commit()

        # Merge data
        merge_table(cursor, table, "srcdb")
        conn.commit()

    conn.execute("DETACH DATABASE srcdb")
    conn.close()

def main():
    
    logging.basicConfig(
        filename=snakemake.log[0],
        level=logging.INFO,
        format="%(asctime)s %(levelname)s:%(message)s"
    )
    input_databases=list(snakemake.input) 
    out_database=snakemake.output.database

    tables_to_merge = ['busco_sequences', 'full_table']

    # Remove output DB if it exists
    if os.path.exists(out_database):
        os.remove(out_database)

    # Initialize the output DB with the schema and data from the first input DB
    logging.info(f"Initializing output database with {input_databases[0]}")
    #print(f"Initializing output database with {input_databases[0]}")
    merge_databases(out_database, input_databases[0], tables_to_merge)

    # Merge the rest of the input databases
    for db_path in input_databases[1:]:
        logging.info(f"Merging database {db_path}")
        #print(f"Merging database {db_path}")
        merge_databases(out_database, db_path, tables_to_merge)

if __name__ == '__main__':
    main()