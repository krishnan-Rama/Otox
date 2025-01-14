#!/usr/bin/env python3

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import argparse
import subprocess
import os
import matplotlib.patches as patches
from collections import defaultdict
from Bio import Phylo
import pandas as pd
import seaborn as sns
import numpy as np

# ------------------- Domain Annotation Parser ------------------- #
def parse_domain_annotation(domain_file):
    """
    Parse hmmsearch domain table output and return structured annotation:
    {sequence: [(domain, start, end, accession, e-value)]}.
    """
    domain_annotations = defaultdict(list)
    with open(domain_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split()
            if len(cols) < 20:
                continue

            seq_name = cols[0]
            domain_name = cols[3]
            accession = cols[4]
            start = int(cols[17])
            end = int(cols[18])
            e_value = float(cols[6])

            domain_annotations[seq_name].append((domain_name, start, end, accession, e_value))

    return domain_annotations

# ------------------- Artistic Domain Visualization ------------------- #
def visualize_domains(domain_annotations, output_image):
    """
    Create a proportional and artistic domain visualization for ortholog sequences.
    """
    sequences = list(domain_annotations.keys())
    max_length = max([max([end for _, _, end, _, _ in domains]) for domains in domain_annotations.values()])

    fig, ax = plt.subplots(figsize=(12, 0.6 * len(sequences)))
    y_position = 0

    domain_colors = {}
    available_colors = plt.cm.tab20.colors
    color_index = 0

    for seq in sorted(sequences):
        domains = sorted(domain_annotations[seq], key=lambda x: x[1])
        y_position += 1
        ax.text(-5, y_position, seq, va='center', ha='right', fontsize=8)

        # Draw sequence line
        ax.add_patch(
            patches.Rectangle((0, y_position - 0.1), max_length, 0.2, facecolor="lightgrey", edgecolor=None)
        )
        # Track motif lanes
        lanes = []

        for domain_name, start, end, accession, _ in domains:
            lane_assigned = False
            for lane in lanes:
                if all(not (existing_start < end and existing_end > start) for existing_start, existing_end in lane):
                    lane.append((start, end))
                    lane_assigned = True
                    break
            if not lane_assigned:
                lanes.append([(start, end)])

            # Assign colors
            if domain_name not in domain_colors:
                domain_colors[domain_name] = (available_colors[color_index % len(available_colors)], accession)
                color_index += 1

            # Draw motif blocks
            for lane_index, lane in enumerate(lanes):
                if (start, end) in lane:
                    y_offset = y_position + lane_index * 0.2  # Keep stacking without gaps
                    ax.add_patch(
                        patches.Rectangle(
                            (start, y_offset - 0.1),
                            end - start,
                            0.2,
                            facecolor=domain_colors[domain_name][0],
                            edgecolor="black",
                            lw=0.5
                        )
                    )
                    # Position residue numbers more legibly above the blocks
                    ax.text(
                        start + (end - start) / 2,
                        y_offset + 0,  # Raise the text for better visibility
                        f"{start}-{end}",
                        fontsize=6,
                        ha="center",
                        color="black"
                    )
                    break

    ax.set_xlim(0, max_length)
    ax.set_ylim(0, y_position + len(lanes) * 0.2)
    ax.set_yticks([])
    ax.set_xlabel("Residue Position", fontsize=10)
    ax.set_title("Domain Annotation on Ortholog Sequences", fontsize=12)

    # Add legend with domain accession numbers
    if domain_colors:
        handles = [
            patches.Patch(color=color, label=f"{domain} ({accession})")
            for domain, (color, accession) in domain_colors.items()
        ]
        ax.legend(
            handles=handles,
            loc='upper right',
            bbox_to_anchor=(1, 0.5),
            fontsize=8,
            title="Domains",
            frameon=True
        )

    # Adjust layout to ensure legend is visible
    plt.tight_layout(rect=[0, 0, 0.85, 1])  # Leave space for legend on the right

    plt.savefig(f"{output_image}.png", dpi=300)
    plt.savefig(f"{output_image}.pdf", dpi=300)
    plt.close()
    print(f"Domain visualization saved to {output_image}.png and {output_image}.pdf")


# ------------------- Domain Conservation Heatmap ------------------- #
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


def visualize_domain_conservation(domain_file, output_image):
    """
    Generate a heatmap for domain conservation levels across ortholog sequences.
    """
    domain_data = []
    with open(domain_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split()
            if len(cols) < 20:
                continue

            seq_name = cols[0]
            domain_name = cols[3]
            accession = cols[4]  # Extract accession number
            c_evalue = float(cols[6])

            domain_data.append([seq_name, domain_name, accession, c_evalue])

    # Convert to DataFrame
    df = pd.DataFrame(domain_data, columns=["Sequence", "Domain", "Accession", "c-Evalue"])

    # Add combined domain label with accession
    df["Domain_Label"] = df["Domain"] + " (" + df["Accession"] + ")"

    # Aggregate to handle duplicates by taking the minimum c-Evalue
    df = df.groupby(["Sequence", "Domain_Label"], as_index=False).agg({"c-Evalue": "min"})

    # Compute natural log-transformed c-Evalues
    df["Log(c-Evalue)"] = -np.log10(df["c-Evalue"])

    # Pivot table for heatmap
    heatmap_data = df.pivot(index="Sequence", columns="Domain_Label", values="Log(c-Evalue)")

    # Define colormap
    cmap = sns.color_palette("viridis", as_cmap=True)  # Blue to Green gradient

    # Create heatmap
    plt.figure(figsize=(16, 12))
    sns.heatmap(
        heatmap_data,
        cmap=cmap,
        cbar_kws={"label": "Conservation (-log10 c-Evalue)"},
        linewidths=0.5,
        linecolor="gray",
        mask=heatmap_data.isnull(),
    )
    plt.title("Domain Conservation Across Ortholog Sequences", fontsize=16, weight="bold")
    plt.xlabel("Domains (with Accession Numbers)", fontsize=14, weight="bold")
    plt.ylabel("Ortholog Sequences", fontsize=14, weight="bold")
    plt.tight_layout()

    # Save the heatmap
    plt.savefig(f"{output_image}.png", dpi=300)
    plt.savefig(f"{output_image}.pdf", dpi=300)
    plt.close()
    print(f"Heatmap saved to {output_image}.png and {output_image}.pdf")

# ------------------- Sequence Alignment ------------------- #
def align_sequences(input_file, output_file):
    """
    Align sequences using MAFFT and save the result to a file.
    """
    print(f"Running MAFFT to align sequences from {input_file}...")
    cmd = ["mafft", "--auto", "--quiet", input_file]
    with open(output_file, "w") as outfile:
        subprocess.run(cmd, stdout=outfile, check=True)
    print(f"Aligned sequences saved to {output_file}")

# ------------------- Phylogenetic Analysis ------------------- #
def construct_phylogenetic_tree(aligned_file, output_tree):
    """
    Construct a phylogenetic tree using FastTree and save it to a file.
    """
    print(f"Constructing phylogenetic tree from {aligned_file}...")
    cmd = ["fasttree", aligned_file]
    with open(output_tree, "w") as outfile:
        subprocess.run(cmd, stdout=outfile, check=True)
    print(f"Phylogenetic tree saved to {output_tree}")

def visualize_phylogenetic_tree(tree_file, output_dir):
    """
    Visualize rooted and unrooted phylogenetic trees dynamically with high-quality layouts.
    """
    try:
        print(f"Reading tree from {tree_file}")
        tree = Phylo.read(tree_file, "newick")

        # Rooted Tree Visualization
        rooted_output = os.path.join(output_dir, "rooted_tree.png")
        plt.figure(figsize=(10, 8))
        Phylo.draw(tree, do_show=False)
        plt.savefig(rooted_output, dpi=300)
        plt.close()
        print(f"Rooted tree saved to {rooted_output}")

        # Unrooted Tree Visualization
        unrooted_output = os.path.join(output_dir, "unrooted_tree.png")
        plt.figure(figsize=(10, 8))
        Phylo.draw(tree, do_show=False)
        plt.savefig(unrooted_output, dpi=300)
        plt.close()
        print(f"Unrooted tree saved to {unrooted_output}")

    except Exception as e:
        print(f"Error during tree visualization: {e}")

# ------------------- Run Domain Annotation ------------------- #
def run_hmmsearch(input_file, output_file):
    """Run HMMER hmmsearch on ungapped sequences."""
    pfam_db = "/opt/conda/envs/bioinformatics/share/Pfam-A.hmm"
    print(f"Running hmmsearch on {input_file}...")
    cmd = ["hmmsearch", "--domtblout", output_file, pfam_db, input_file]
    subprocess.run(cmd, check=True)

def ungap_fasta(input_file, output_file):
    """Remove gaps from aligned sequences for hmmsearch."""
    with open(input_file) as infile, open(output_file, "w") as outfile:
        seq = ""
        header = None
        for line in infile:
            if line.startswith(">"):
                if header:
                    outfile.write(header + "\n" + seq + "\n")
                header = line.strip()
                seq = ""
            else:
                seq += line.strip().replace("-", "")
        if header:
            outfile.write(header + "\n" + seq + "\n")

# ------------------- Main Function ------------------- #
def main():
    parser = argparse.ArgumentParser(description="Domain Annotation and Visualization Pipeline")
    parser.add_argument("--input", required=True, help="Input FASTA file (unaligned orthologs)")
    parser.add_argument("--output", default="results", help="Output directory")
    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    # Filenames
    unaligned_file = args.input
    aligned_file = os.path.join(args.output, "aligned_sequences.fa")
    ungapped_file = os.path.join(args.output, "ungapped_sequences.fa")
    domain_outfile = os.path.join(args.output, "domain_annotation.domtblout")
    domain_plot = os.path.join(args.output, "domain_visualization")
    domain_heatmap = os.path.join(args.output, "domain_conservation_heatmap")
    tree_file = os.path.join(args.output, "phylogenetic_tree.nwk")

    # 1. Perform sequence alignment
    print("Aligning sequences...")
    align_sequences(unaligned_file, aligned_file)

    # 2. Construct phylogenetic tree
    print("Constructing phylogenetic tree...")
    construct_phylogenetic_tree(aligned_file, tree_file)

    # 3. Ungap sequences for HMMER
    print("Removing gaps from sequences for HMMER...")
    ungap_fasta(unaligned_file, ungapped_file)

    # 4. Run HMMER hmmsearch
    print("Running domain annotation...")
    run_hmmsearch(ungapped_file, domain_outfile)

    # 5. Parse domain annotation
    print("Parsing domain annotations...")
    domain_annotations = parse_domain_annotation(domain_outfile)

    # 6. Visualize domain annotations
    print("Creating domain annotation visualization...")
    visualize_domains(domain_annotations, domain_plot)

    # 7. Create domain conservation heatmap
    print("Creating domain conservation heatmap...")
    visualize_domain_conservation(domain_outfile, domain_heatmap)

    print("Pipeline completed. Results are saved in:", args.output)

if __name__ == "__main__":
    main()
