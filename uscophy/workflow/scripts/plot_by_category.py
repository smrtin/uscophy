import click
import pandas as pd
import sqlite3
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import math
import sys


@click.command()
@click.option('--sample_file', type=click.Path(exists=True), required=True,
              help='Sample-category file with Sample as first column and one or more category columns')
@click.option('--sqlite_file', type=click.Path(exists=True), required=True,
              help='SQLite3 database file with BUSCO results')
@click.option('--output', type=click.Path(), default='status_boxplots.pdf', show_default=True,
              help='Output PDF file for plots')
@click.option('--categories', multiple=True,
              help='Specific categories (column headers) to plot; default is all category columns')
def main(sample_file, sqlite_file, output, categories):
    """
    Generates boxplots of unique BUSCO Status counts per sample and category column.
    Supports multiple category columns in the sample file.
    Produces one combined plot for all samples,
    then faceted plots for each specified category column.
    """
    # Read sample-category table
    try:
        sample_df = pd.read_csv(sample_file, sep=None, engine='python', header=0)
    except Exception as e:
        click.secho(f"Error reading sample file: {e}", fg='red')
        sys.exit(1)

    # Validate presence of 'Sample' column
    if 'Sample' not in sample_df.columns:
        click.secho(f"'Sample' column not found in sample file", fg='red')
        sys.exit(1)

    # Extract category columns (all except 'Sample')
    category_cols = [col for col in sample_df.columns if col != 'Sample']

    if not category_cols:
        click.secho(f"No category columns found beyond 'Sample' in sample file", fg='red')
        sys.exit(1)

    # If user specified categories, validate and use only those
    if categories:
        selected_categories = [cat for cat in categories if cat in category_cols]
        missing_cats = set(categories) - set(selected_categories)
        if missing_cats:
            click.secho(f"Warning: Some specified categories not found and will be skipped: {', '.join(missing_cats)}", fg='yellow')
        if not selected_categories:
            click.secho("No valid categories specified to plot - exiting.", fg='red')
            sys.exit(1)
    else:
        selected_categories = category_cols  # Plot all by default

    # Fill missing category values with 'NoCategory' for all category columns to plot
    sample_df[selected_categories] = sample_df[selected_categories].fillna('NoCategory')

    # Load BUSCO results DB
    try:
        conn = sqlite3.connect(sqlite_file)
        db_df = pd.read_sql_query('SELECT Species, Status, Busco_id FROM full_table', conn)
        conn.close()
    except Exception as e:
        click.secho(f"Error reading SQLite database: {e}", fg='red')
        sys.exit(1)

    samples_in_file = set(sample_df['Sample'])
    samples_in_db = set(db_df['Species'])
    missing_samples = sorted(samples_in_db - samples_in_file)
    if missing_samples:
        click.secho("NOTE: Samples in database with no category (not in sample_file):\n" +
                    ', '.join(str(x) for x in missing_samples if x is not None),
                    fg='cyan')

    # Filter db for samples only present in sample file
    db_df = db_df[db_df['Species'].isin(samples_in_file)]

    # Count unique Busco_id per (Species, Status)
    statuses = ['Complete', 'Fragmented', 'Duplicated', 'Missing']
    unique_counts = (
        db_df.groupby(['Species', 'Status'])['Busco_id']
        .nunique()
        .unstack(fill_value=0)
        .reindex(columns=statuses, fill_value=0)
    )

    # Prepare combined plot data (all samples, no category splitting)
    combined = unique_counts.reset_index().rename(columns={'Species': 'Sample'})
    combined_melted = combined.melt(
        id_vars=['Sample'],
        value_vars=statuses,
        var_name='Status',
        value_name='Count'
    )
    ymax = combined_melted['Count'].max() * 1.05 if combined_melted['Count'].max() > 0 else 1

    # Start PDF output
    with PdfPages(output) as pdf:

        # --- Combined plot for all samples ---
        plt.figure(figsize=(8, 6))
        sns.boxplot(
            x='Status', y='Count', hue='Status',
            data=combined_melted, palette='colorblind', legend=False
        )
        sns.stripplot(
            x='Status', y='Count',
            data=combined_melted, color='black', size=3, jitter=True, dodge=True
        )
        plt.title('Unique BUSCO Status Counts per Sample (All Samples)')
        plt.ylim(0, ymax)
        plt.tight_layout()
        pdf.savefig()
        plt.close()

        # --- Faceted plots by each selected category ---
        for category in selected_categories:
            merged = unique_counts.reset_index().merge(
                sample_df[['Sample', category]], left_on='Species', right_on='Sample', how='left'
            )
            merged[category] = merged[category].fillna('NoCategory')

            melted = merged.melt(
                id_vars=['Sample', category],
                value_vars=statuses,
                var_name='Status',
                value_name='Count'
            ).rename(columns={category: 'Category'})

            # Sort categories for consistent plotting, 'NoCategory' last
            cats = sorted([c for c in melted['Category'].unique() if c != 'NoCategory'])
            if 'NoCategory' in melted['Category'].values:
                cats.append('NoCategory')
            melted['Category'] = pd.Categorical(melted['Category'], categories=cats, ordered=True)

            n_cat = len(cats)
            col_wrap = math.ceil(math.sqrt(n_cat))

            g = sns.catplot(
                x='Status', y='Count', col='Category',
                data=melted, kind='box', hue='Status',
                palette='colorblind', legend=False,
                col_wrap=col_wrap, sharey=False
            )
            for ax, cat in zip(g.axes.flatten(), cats):
                sns.stripplot(
                    x='Status', y='Count',
                    data=melted[melted['Category'] == cat],
                    color='black', size=3, jitter=True, dodge=True, ax=ax
                )
                ax.set_ylim(0, ymax)
            g.fig.subplots_adjust(top=0.85)
            g.fig.suptitle(f'Unique BUSCO Status Counts by Sample Category ({category})')

            pdf.savefig(g.fig)
            plt.close('all')

    click.echo(f"Plots saved to {output}")


if __name__ == '__main__':
    main()
