import sqlite3
import matplotlib.pyplot as plt
import numpy as np
import click
import matplotlib.ticker as ticker

@click.command()
@click.argument("database", type=click.Path(exists=True))
@click.option("--output", default="species_busco_threshold_barplot.png", show_default=True, 
              help="Output PNG file name for the plot")


def plot_species_threshold_lineplot(database, output):
    """
    Plot number of species retained at varying BUSCO thresholds from DATABASE,
    using the fixed table 'full_table', and save the plot to OUTPUT file as a line plot.
    """
    with sqlite3.connect(database) as db_connection:
        cursor = db_connection.cursor()

        table = "full_table"

        # Get total number of unique BUSCO ids
        cursor.execute(f"SELECT COUNT(DISTINCT Busco_id) FROM {table}")
        total_unique_busco = cursor.fetchone()[0]
        click.echo(f"Total unique BUSCO IDs: {total_unique_busco}")

        # Get all species
        cursor.execute(f"SELECT DISTINCT Species FROM {table}")
        all_species = [row[0] for row in cursor.fetchall()]
        click.echo(f"Found {len(all_species)} species in the table")

        xs = []
        complete_ys = []
        complete_frag_ys = []
        complete_frag_dup_ys = []
        
        for p in range(5, 101, 5):
            min_genes = int(np.ceil(total_unique_busco * (p / 100)))
            
            # Only "Complete"
            pass_count_complete = 0
            for species in all_species:
                cursor.execute(
                    f"SELECT COUNT(DISTINCT Busco_id) FROM {table} "
                    "WHERE Species=? AND Status='Complete'",
                    (species,)
                )
                n = cursor.fetchone()[0]
                if n >= min_genes:
                    pass_count_complete += 1
            
            # "Complete" or "Fragmented"
            pass_count_frag = 0
            for species in all_species:
                cursor.execute(
                    f"SELECT COUNT(DISTINCT Busco_id) FROM {table} "
                    "WHERE Species=? AND (Status='Complete' OR Status='Fragmented')",
                    (species,)
                )
                n = cursor.fetchone()[0]
                if n >= min_genes:
                    pass_count_frag += 1
            
            # "Complete" or "Fragmented"
            pass_count_dup = 0
            for species in all_species:
                cursor.execute(
                    f"SELECT COUNT(DISTINCT Busco_id) FROM {table} "
                    "WHERE Species=? AND (Status='Complete' OR Status='Fragmented' OR Status='Duplicated')",
                    (species,)
                )
                n = cursor.fetchone()[0]
                if n >= min_genes:
                    pass_count_dup += 1

            xs.append(p)
            complete_ys.append(pass_count_complete)
            complete_frag_ys.append(pass_count_frag)
            complete_frag_dup_ys.append(pass_count_dup)
            click.echo(
                f"Threshold {p}%: min_genes={min_genes}, Complete={pass_count_complete}, Complete+Frag={pass_count_frag}, Complete+Frag+Dupl={pass_count_dup}"
            )

    fig, ax = plt.subplots()
    ax.plot(xs, complete_ys, marker='o', linestyle='-', alpha=0.6, label="Complete only")
    ax.plot(xs, complete_frag_ys, marker='s', linestyle='-', alpha=0.6, label="Complete or Fragmented")
    ax.plot(xs, complete_frag_dup_ys, marker='s', linestyle='-', alpha=0.6, label="Complete or Fragmented or best Duplicate")

    ax.set_xlabel("Minimum BUSCOs present (%)")
    ax.set_ylabel("Number of species retained")
    ax.set_title("Species retained at varying BUSCO presence thresholds")
    ax.legend()
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))  # Set ticks every 10%
    ax.set_xlim(0, 105)  # A bit past 100 to nicely include last tick
    plt.tight_layout()
    plt.savefig(output)
    click.echo(f"Plot saved to {output}")

if __name__ == "__main__":
    plot_species_threshold_lineplot()
