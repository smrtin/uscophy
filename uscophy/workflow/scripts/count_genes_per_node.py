import pandas as pd
from ete3 import Tree
import click

def standardize(name):
    return name.strip().replace(' ', '_')

def read_table(filename):
    with open(filename, 'r', encoding='utf-8') as f:
        sample = f.read(4096)
        import csv
        sniffer = csv.Sniffer()
        dialect = sniffer.sniff(sample, delimiters=[',', '\t'])
        delimiter = dialect.delimiter
    df = pd.read_csv(filename, delimiter=delimiter, header=0, index_col=0)
    df.columns = [standardize(col) for col in df.columns]
    return df

def read_tree(newick_file):
    t = Tree(newick_file, format=1)
    for leaf in t.get_leaves():
        leaf.name = standardize(leaf.name)
    return t

def root_tree(tree, root_species):
    if root_species:
        std_root_species = standardize(root_species)
        leaf = tree.search_nodes(name=std_root_species)
        if not leaf:
            raise ValueError(f"Root species '{root_species}' not found in tree tips.")
        tree.set_outgroup(leaf[0])

def match_tree_and_traits(tree, gene_df):
    tree_tips = [leaf.name for leaf in tree.get_leaves()]
    gene_df = gene_df[[tip for tip in tree_tips if tip in gene_df.columns]]
    gene_df = gene_df.loc[:, ~gene_df.columns.duplicated()]
    return gene_df

def fitch_pass(tree, gene_states):
    # Assign states to leaves
    for leaf in tree.iter_leaves():
        state = gene_states.get(leaf.name, None)
        if state is None or pd.isnull(state):
            leaf.add_feature("fitch_set", set(["0", "1"]))  # ambiguous if missing
        else:
            leaf.add_feature("fitch_set", set([str(int(state))]))
    # Upward pass
    for node in tree.traverse("postorder"):
        if not node.is_leaf():
            child_sets = [child.fitch_set for child in node.get_children()]
            inter = set.intersection(*child_sets)
            if inter:
                node.add_feature("fitch_set", inter)
            else:
                node.add_feature("fitch_set", set.union(*child_sets))
    # Downward pass (assign state: pick "1" if possible, else "0")
    def assign_state(node, parent_state=None):
        if parent_state and parent_state in node.fitch_set:
            node.add_feature("fitch_state", parent_state)
        else:
            node.add_feature("fitch_state", "1" if "1" in node.fitch_set else "0")
        for child in node.get_children():
            assign_state(child, node.fitch_state)
    assign_state(tree)

def annotate_tree_with_gene_counts(tree, gene_df):
    # For each gene, reconstruct presence/absence at all nodes
    for gene in gene_df.index:
        char_states = {col: gene_df.loc[gene, col] for col in gene_df.columns}
        fitch_pass(tree, char_states)
        for node in tree.traverse():
            if not hasattr(node, "gene_states"):
                node.gene_states = {}
            node.gene_states[gene] = node.fitch_state
    # Count number of genes present at each node and annotate
    for idx, node in enumerate(tree.traverse()):
        count = sum(1 for gene in gene_df.index if node.gene_states.get(gene, "0") == "1")
        if node.is_leaf():
            base_name = node.name.split("__genes_")[0]
            node.name = f"{base_name}__genes_{count}"
        else:
            node.name = f"Node{idx}__genes_{count}"

@click.command()
@click.option('--table', '-t', required=True, help='Gene presence/absence table CSV/TSV file')
@click.option('--tree', '-n', required=True, help='Newick tree file')
@click.option('--output', '-o', default='annotated_tree.nwk', show_default=True, help='Output filename for the annotated Newick tree')
@click.option('--root-species', '-r', default=None, help='Species name to use as outgroup for rooting the tree')


def main(table, tree, output, root_species):
    gene_df = read_table(table)
    t = read_tree(tree)
    root_tree(t, root_species)
    gene_df = match_tree_and_traits(t, gene_df)
    annotate_tree_with_gene_counts(t, gene_df)
    t.write(outfile=output, format=1)  # format=1: include internal node names
    print(f"Annotated tree written to: {output}")

if __name__ == "__main__":
    main()
