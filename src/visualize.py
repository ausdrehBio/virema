import argparse
from pathlib import Path

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


##########  Deletion (Forward)  ##########
# This occurs when the polymerase skips forward along the template,
# leaving out a section of the genome.
# The Mechanism:
# The polymerase copies up to a certain point (e.g., nucleotide 300), "falls off" or jumps forward,
# and lands further downstream (e.g., nucleotide 1100), resuming copying from there.
# The Result:
# The resulting RNA is shorter than the original because the segment between 300 and 1100 is missing.
# How to spot it:
# The Start coordinate is smaller than the Stop coordinate ($Start < Stop$).
# Example from your data:
# 303_to_1102 (Start 303 is less than Stop 1102).

##########  Duplication (Reverse)   ##########
# This occurs when the polymerase jumps backward to a section it has already copied.
# The Mechanism:
# The polymerase copies up to a point (e.g., nucleotide 1200), but then jumps back to an earlier point
# (e.g., nucleotide 1100) and copies that same section again.
# The Result:
# The resulting RNA is longer than normal because the segment between 1100 and 1200 is repeated (duplicated).
# How to spot it:
# The Start coordinate is larger than the Stop coordinate ($Start > Stop$).
# Example from your data:
# 1153_to_1132 (Start 1153 is greater than Stop 1132).


def _safe_filename(text):
    return "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in text)


def parse_virema_bed(filepath: Path):
    """Parses the ViReMa BED file to extract coordinates and classifications."""
    if not filepath.exists():
        print(f"Error: File not found at {filepath}")
        return pd.DataFrame()

    data = []

    with filepath.open("r") as f:
        for line in f:
            line = line.strip()
            # Skip header lines (BED files often start with 'track')
            if line.startswith("track") or not line:
                continue

            # ViReMa BED Format usually:
            # Chrom | Start | Stop | Type | Count | Strand | ...
            cols = line.split("\t")

            try:
                chrom = cols[0]
                start = int(cols[1])
                stop = int(cols[2])
                event_type = cols[3]
                count = int(cols[4])

                data.append(
                    {
                        "Chrom": chrom,
                        "Start": start,
                        "Stop": stop,
                        "Type": event_type,
                        "Count": count,
                        "Jump_Size": abs(stop - start),
                    }
                )
            except (ValueError, IndexError):
                continue

    return pd.DataFrame(data)


def plot_arc_diagram(df, chrom_name, output_dir, genome_length=None, show=False):
    """Plots arcs colored by the Classification found in the BED file."""
    subset = df[df["Chrom"] == chrom_name]

    if subset.empty:
        print(f"No data found for chromosome: {chrom_name}")
        return

    fig, ax = plt.subplots(figsize=(14, 6))
    max_count = subset["Count"].max()

    print(f"Plotting {len(subset)} events for {chrom_name}...")

    color_map = {
        "Deletion": "red",
        "Duplication": "blue",
        "Back-Splice": "purple",
        "Insertion": "green",
        "Splice": "orange",
    }

    for _, row in subset.iterrows():
        start, stop = row["Start"], row["Stop"]
        center = (start + stop) / 2
        width = abs(stop - start)
        height = width / 2

        if max_count <= 1:
            alpha = 1.0
        else:
            alpha = min(1.0, 0.1 + 0.9 * (np.log1p(row["Count"]) / np.log1p(max_count)))

        event_type = row["Type"]
        color = "gray"
        for key, c in color_map.items():
            if key in event_type:
                color = c
                break

        arc = patches.Arc(
            (center, 0),
            width,
            height * 2,
            theta1=0,
            theta2=180,
            edgecolor=color,
            alpha=alpha,
            linewidth=1.5,
        )
        ax.add_patch(arc)

    limit = subset[["Start", "Stop"]].max().max()
    if genome_length:
        limit = max(limit, genome_length)
    ax.set_xlim(0, limit)
    ax.set_ylim(0, max(1, subset["Jump_Size"].max() / 1.5))
    ax.set_xlabel("Genome Position (nt)")
    ax.set_ylabel("Jump Distance")
    ax.set_title(f"Recombination Arc Diagram (BED Data)\n{chrom_name}")

    legend_elements = [Line2D([0], [0], color=c, lw=2, label=t) for t, c in color_map.items()]
    ax.legend(handles=legend_elements, loc="upper right", title="Event Type")

    plt.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    filename = f"arc_{_safe_filename(chrom_name)}.png"
    fig.savefig(output_dir / filename, dpi=300)
    if show:
        plt.show()
    plt.close(fig)


def plot_top_events_bed(df, chrom_name, output_dir, top_n=15, show=False):
    """Bar chart of top events using BED classifications."""
    subset = df[df["Chrom"] == chrom_name].copy()
    if subset.empty:
        return

    subset["Label"] = (
        subset["Start"].astype(str)
        + "->"
        + subset["Stop"].astype(str)
        + "\n("
        + subset["Type"]
        + ")"
    )

    top_events = subset.nlargest(top_n, "Count")

    fig, ax = plt.subplots(figsize=(12, 6))
    colors = top_events["Type"].apply(
        lambda x: "red" if "Deletion" in x else ("blue" if "Duplication" in x else "gray")
    )

    bars = ax.bar(top_events["Label"], top_events["Count"], color=colors)

    ax.set_ylabel("Read Count")
    ax.set_title(f"Top {top_n} Events by Classification\n{chrom_name}")
    plt.xticks(rotation=45, ha="right")

    for bar in bars:
        height = bar.get_height()
        ax.annotate(
            f"{height}",
            xy=(bar.get_x() + bar.get_width() / 2, height),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
        )

    plt.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    filename = f"top_events_{_safe_filename(chrom_name)}.png"
    fig.savefig(output_dir / filename, dpi=300)
    if show:
        plt.show()
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(description="Visualize ViReMa BED results.")
    parser.add_argument("bed_file", help="Path to Virus_Recombination_Results.bed")
    parser.add_argument("--output-dir", required=True, help="Directory to write plots.")
    parser.add_argument("--chrom", action="append", help="Chrom/segment to plot (repeatable).")
    parser.add_argument("--genome-length", type=int, help="Optional genome length for axis limits.")
    parser.add_argument("--top-n", type=int, default=15, help="Number of top events to plot.")
    parser.add_argument("--show", action="store_true", help="Display plots interactively.")
    return parser.parse_args()


def main():
    args = parse_args()
    bed_path = Path(args.bed_file)
    output_dir = Path(args.output_dir)

    print(f"Reading BED file: {bed_path}...")
    df = parse_virema_bed(bed_path)
    if df.empty:
        print("Error: DataFrame is empty. Check if the BED file exists and has data.")
        return

    chroms = args.chrom or sorted(df["Chrom"].unique())
    print(f"Found {len(chroms)} chromosome(s) to visualize.")
    for chrom in chroms:
        plot_arc_diagram(
            df,
            chrom,
            output_dir,
            genome_length=args.genome_length,
            show=args.show,
        )
        plot_top_events_bed(
            df,
            chrom,
            output_dir,
            top_n=args.top_n,
            show=args.show,
        )


if __name__ == "__main__":
    main()
