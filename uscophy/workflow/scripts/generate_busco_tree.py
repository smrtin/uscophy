#!/usr/bin/env python3

import click
import subprocess
import os
import glob
import csv
import re
import pathlib

def parse_busco_summary(summary_file):
    """
    Parse a BUSCO summary file and return counts for each category and lineage.
    Args:
        summary_file (str): Path to BUSCO summary text file.
    Returns:
        values (list): BUSCO counts [complete_single, complete_duplicated, fragmented, missing].
        lineage (str or None): BUSCO lineage name, detected from filename or summary content.
    """
    categories = [
        'Complete and single-copy BUSCOs',
        'Complete and duplicated BUSCOs',
        'Fragmented BUSCOs',
        'Missing BUSCOs'
    ]
    values = [0, 0, 0, 0]
    lineage = None

    # Attempt to extract lineage from filename pattern
    basename = os.path.basename(summary_file)
    match = re.match(r"short_summary\.specific\.([^.]+)\..+\.txt", basename)
    if match:
        lineage = match.group(1)

    # Parse BUSCO counts and lineage info inside summary file
    with open(summary_file) as f:
        for line in f:
            line_strip = line.strip()
            for i, cat in enumerate(categories):
                if cat in line_strip:
                    nums = re.findall(r'\d+', line_strip)
                    if nums:
                        values[i] = int(nums[0])
            if lineage is None and 'lineage dataset is:' in line.lower():
                m = re.search(r'lineage dataset is:?\s*([\w_]+)', line.lower())
                if m:
                    lineage = m.group(1)
    return values, lineage

@click.command()
@click.option('--newick', '-n', required=True, type=click.Path(exists=True),
              help='Path to the Newick tree file.')
@click.option('--busco-folder', '-b', required=True, type=click.Path(exists=True, file_okay=False),
              help='Directory containing BUSCO summary files.')
@click.option('--output', '-o', default="busco_tree_plot", type=click.Path(),
              help='Basename for output plot (with optional extension, default pdf).')
@click.option('--outgroup', '-g', default=None,
              help='Taxon name to use as outgroup for rooting the tree.')
@click.option('--annotation', '-a', required=False, type=click.Path(exists=True),
              help='[Optional] CSV file with first column "Sample" and other columns as annotation categories.')
def run_busco_plot(newick, busco_folder, output, outgroup, annotation):
    """
    Main function to parse BUSCO summaries, process annotations (if provided),
    generate a multi-page PDF plot with one page per annotation category (or single plot without annotation).
    """

    # Locate BUSCO summary files (.txt and .summary)
    summary_files = sorted(
        glob.glob(os.path.join(busco_folder, "*.txt")) +
        glob.glob(os.path.join(busco_folder, "*.summary"))
    )
    if not summary_files:
        raise FileNotFoundError(f"No BUSCO summary files found in {busco_folder}")

    # Parse BUSCO summaries and collect lineage names alongside counts
    all_lineages = set()
    summary_rows = []
    for summary_file in summary_files:
        filename = os.path.basename(summary_file)
        name_parts = filename.split('.')
        # Extract sample name depending on BUSCO naming conventions
        if len(name_parts) > 3 and name_parts[0].startswith("short_summary"):
            taxa = name_parts[-2]
        else:
            taxa = os.path.splitext(filename)[0]
        values, lineage = parse_busco_summary(summary_file)
        if lineage:
            all_lineages.add(lineage)
        summary_rows.append([taxa] + values)

    # Compose lineage title for the legend, handle multiple lineages gracefully
    if len(all_lineages) <= 1:
        lineage_title = next(iter(all_lineages), '')
    else:
        lineage_title = ", ".join(sorted(all_lineages))

    # Write the combined BUSCO summary for R plotting
    csv_filename = "busco_summary.csv"
    with open(csv_filename, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["taxa", "Complete_single_copy", "Complete_duplicated", "Fragmented", "Missing"])
        writer.writerows(summary_rows)

    # Prepare output filename and extension
    output_path = pathlib.Path(output)
    if output_path.suffix:
        base_output = str(output_path.with_suffix(''))
        file_ext = output_path.suffix.lstrip('.')
    else:
        base_output = str(output)
        file_ext = "pdf"
    plot_filename = f"{base_output}.{file_ext}"

    # Process annotation file if provided; extract categories
    if annotation:
        with open(annotation) as annfile:
            reader = csv.DictReader(annfile)
            headers = reader.fieldnames
            if not headers or headers[0] != "Sample":
                raise Exception("Annotation file must have 'Sample' as the first column.")
            categories = headers[1:]
            if not categories:
                raise Exception("Annotation file must have at least one annotation category (columns after 'Sample').")
        categories_r = ", ".join(f'"{cat}"' for cat in categories)
    else:
        categories = []
        categories_r = '""'  # Placeholder for R code when no annotation

    # Generate R script dynamically to handle both annotated and non-annotated workflow,
    # including dynamic height based on number of taxa
    r_script = f'''
library(ggtree)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(ape)
library(patchwork)
library(grid)
library(viridis)

tree_file <- "{newick}"
busco_csv <- "{csv_filename}"
output_file <- "{plot_filename}"
outgroup <- {"NULL" if outgroup is None else f'"{outgroup}"'}
categories <- c({categories_r})
annotation_file <- {('"' + annotation + '"') if annotation else "NULL"}

tree <- read.tree(tree_file)
if (!is.null(outgroup) && outgroup %in% tree$tip.label) {{
  tree <- root(tree, outgroup=outgroup, resolve.root=TRUE)
}}

# Set plot height dynamically according to number of taxa
n_taxa <- length(tree$tip.label)
height <- min(0.2 * n_taxa + 3, 50)
width <- 18  # Fixed width to match previous default

busco <- read_csv(busco_csv, show_col_types = FALSE)
busco$taxa <- as.character(busco$taxa)

if (!is.null(annotation_file)) {{
  ann <- read.csv(annotation_file, stringsAsFactors=FALSE)
  colnames(ann)[1] <- "Sample"
}} else {{
  ann <- NULL
}}

pdf(output_file, width=width, height=height, onefile=TRUE)

if (is.null(ann)) {{
  ggtr <- ggtree(tree) + geom_tiplab(align=TRUE)

  tree_data <- ggtr$data %>% filter(isTip)
  busco_long <- pivot_longer(
    busco,
    cols = c(Complete_single_copy, Complete_duplicated, Fragmented, Missing),
    names_to = "category",
    values_to = "count"
  )
  busco_long$category <- factor(busco_long$category,
                               levels = c("Complete_single_copy", "Complete_duplicated", "Fragmented", "Missing"))

  busco_long <- left_join(busco_long, tree_data %>% select(label, y), by = c("taxa" = "label"))

  y_limits <- range(tree_data$y) + c(-0.5, 0.5)
  max_x <- max(ggtr$data$x, na.rm = TRUE)
  tip_space <- 0.7 * max_x

  ggtr_aligned <- ggtr +
    scale_y_continuous(limits=y_limits, expand=c(0,0)) +
    xlim(0, max_x + tip_space) +
    theme(
      plot.margin=margin(5, 1, 5, 5),
      axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
    )

  barplot <- ggplot(busco_long, aes(y = y, x = count, fill = category)) +
    geom_col(width=0.8, position=position_stack(reverse=TRUE), orientation="y") +
    scale_fill_manual(values=c(
      "Complete_single_copy"="#009E73",
      "Complete_duplicated"="#E69F00",
      "Fragmented"="#56B4E9",
      "Missing"="#333333"
    ),
    breaks=c("Complete_single_copy", "Complete_duplicated", "Fragmented", "Missing")) +
    guides(fill=guide_legend(title="{lineage_title}")) +
    scale_y_continuous(limits=y_limits, breaks=tree_data$y, labels=tree_data$label, expand=c(0,0)) +
    theme_minimal() +
    theme(
      axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      panel.grid.major.y=element_blank(),
      legend.position="right",
      legend.direction="vertical",
      legend.key.height=unit(0.8, "cm"),
      legend.key.width=unit(1.5, "cm"),
      plot.margin=margin(5,5,5,5)
    ) +
    ylab(NULL) + xlab("BUSCO count")

  final_plot <- ggtr_aligned + barplot + plot_layout(ncol=2, widths=c(2,1))
  grid.newpage()
  grid.text("BUSCO summary", x=0.5, y=0.97, gp=gpar(fontsize=24, fontface="bold"))
  print(final_plot, newpage=FALSE)

}} else {{
  for (category in categories) {{
    ann$Category <- ann[[category]]
    ann$Category <- as.factor(ann$Category)
    category_levels <- levels(ann$Category)
    n_cat <- length(category_levels)
    cat_colors <- setNames(viridis(n_cat), category_levels)
    legend_title_cat <- gsub("_", " ", category)

    ggtr <- ggtree(tree) %<+% ann +
      geom_tippoint(aes(color=Category), size=4, show.legend=TRUE, na.rm=TRUE) +
      geom_tiplab(color="black", align=TRUE, size=4, show.legend=FALSE) +
      scale_color_manual(values=cat_colors, name=legend_title_cat)

    tree_data <- ggtr$data %>% filter(isTip)
    busco_long <- pivot_longer(
      busco,
      cols = c(Complete_single_copy, Complete_duplicated, Fragmented, Missing),
      names_to = "category",
      values_to = "count"
    )
    busco_long$category <- factor(busco_long$category,
                                 levels=c("Complete_single_copy", "Complete_duplicated", "Fragmented", "Missing"))
    busco_long <- left_join(busco_long, tree_data %>% select(label, y), by=c("taxa"="label"))

    y_limits <- range(tree_data$y) + c(-0.5, 0.5)
    max_x <- max(ggtr$data$x, na.rm=TRUE)
    tip_space <- 0.7 * max_x

    ggtr_aligned <- ggtr +
      scale_y_continuous(limits=y_limits, expand=c(0,0)) +
      xlim(0, max_x + tip_space) +
      theme(
        plot.margin=margin(5,1,5,5),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
      )

    barplot <- ggplot(busco_long, aes(y=y, x=count, fill=category)) +
      geom_col(width=0.8, position=position_stack(reverse=TRUE), orientation="y") +
      scale_fill_manual(values=c(
        "Complete_single_copy"="#009E73",
        "Complete_duplicated"="#E69F00",
        "Fragmented"="#56B4E9",
        "Missing"="#333333"
      ),
      breaks=c("Complete_single_copy", "Complete_duplicated", "Fragmented", "Missing")) +
      guides(fill=guide_legend(title="{lineage_title}")) +
      scale_y_continuous(limits=y_limits, breaks=tree_data$y, labels=tree_data$label, expand=c(0,0)) +
      theme_minimal() +
      theme(
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y=element_blank(),
        legend.position="right",
        legend.direction="vertical",
        legend.key.height=unit(0.8, "cm"),
        legend.key.width=unit(1.5, "cm"),
        plot.margin=margin(5,5,5,5)
      ) +
      ylab(NULL) + xlab("BUSCO count")

    final_plot <- ggtr_aligned + barplot + plot_layout(ncol=2, widths=c(2,1))

    category_title <- ifelse(is.null(category) || category == "", "BUSCO summary", category)
    grid.newpage()
    grid.text(category_title, x=0.5, y=0.97, gp=gpar(fontsize=24, fontface="bold"))
    print(final_plot, newpage=FALSE)
  }}
}}

dev.off()
cat("Plots saved to", output_file, "\\n")
'''

    r_script_file = "busco_plot_script.R"
    with open(r_script_file, "w") as rfile:
        rfile.write(r_script)

    try:
        subprocess.run(
            ["Rscript", r_script_file],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True
        )
    except subprocess.CalledProcessError as e:
        print("Error: R script failed.")
        print(e.stderr)
    finally:
        if os.path.exists(r_script_file):
            os.remove(r_script_file)

if __name__ == "__main__":
    run_busco_plot()
