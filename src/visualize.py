import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np


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
# 303_to_1102 (Start 303 is less than Stop 1102).2. Back-splice

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




# --- Configuration ---
INPUT_FILE = "data/_CompiledVirus_Recombination_Results.txt"
# Set this to the specific library you want to visualize (copy exact name from file)
TARGET_LIBRARY = "NC_004146.1_FHV_RNA1.seq_to_NC_004146.1_FHV_RNA1.seq"
GENOME_LENGTH = 3107 # Length of FHV RNA1 (approx standard, adjust if needed)
OUTPUT_DIR = "output"

def _safe_filename(text):
    return "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in text)

def parse_virema_results(filepath):
    """Parses the specific ViReMa compiled format into a Pandas DataFrame."""
    data = []
    current_library = None
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            
            if line.startswith("@NewLibrary:"):
                current_library = line.split(": ")[1].strip()
                continue
            if line.startswith("@EndofLibrary"):
                current_library = None
                continue
                
            # Parse the data lines. Example: 303_to_1102_#_325
            # Lines are tab-separated list of events
            events = line.split('\t')
            for event in events:
                if "_#_" in event:
                    try:
                        # Split 'Start_to_Stop' and 'Count'
                        loc_part, count_part = event.split('_#_')
                        count = int(count_part)
                        
                        # Split 'Start' and 'Stop'
                        start, stop = map(int, loc_part.split('_to_'))
                        
                        data.append({
                            'Library': current_library,
                            'Start': start,
                            'Stop': stop,
                            'Count': count,
                            'Jump_Size': abs(stop - start)
                        })
                    except ValueError:
                        continue
                        
    return pd.DataFrame(data)

def plot_arc_diagram(df, library_name):
    """Plots recombination events as arcs along the genome."""
    subset = df[df['Library'] == library_name]
    
    if subset.empty:
        print(f"No data found for library: {library_name}")
        return

    fig, ax = plt.subplots(figsize=(14, 6))
    
    # Normalize count for opacity (alpha)
    max_count = subset['Count'].max()
    
    print(f"Plotting {len(subset)} unique events for Arc Diagram...")
    
    for _, row in subset.iterrows():
        start, stop = row['Start'], row['Stop']
        center = (start + stop) / 2
        width = abs(stop - start)
        height = width / 2  # Semi-circle height is half width
        
        # Calculate alpha (opacity) based on log scale so low counts are still visible
        alpha = min(1.0, 0.1 + 0.9 * (np.log(row['Count']) / np.log(max_count)))
        
        # Color: Red for deletions (Forward), Blue for Back-splices (Reverse)
        color = 'red' if start < stop else 'blue'
        
        arc = patches.Arc((center, 0), width, height * 2, 
                          theta1=0, theta2=180, edgecolor=color, alpha=alpha, linewidth=1.5)
        ax.add_patch(arc)

    ax.set_xlim(0, max(subset[['Start', 'Stop']].max().max(), GENOME_LENGTH))
    ax.set_ylim(0, subset['Jump_Size'].max() / 1.5)
    ax.set_xlabel("Genome Position (nt)")
    ax.set_ylabel("Arc Height (Jump Distance)")
    ax.set_title(f"Recombination Arc Diagram\n{library_name}")
    ax.text(0.02, 0.95, "Red = Deletion (Forward)\nBlue = Back-splice/Duplication (Reverse)\nOpacity = Frequency", 
            transform=ax.transAxes, fontsize=10, verticalalignment='top')
    
    plt.tight_layout()
    filename = f"arc_{_safe_filename(library_name)}.png"
    fig.savefig(os.path.join(OUTPUT_DIR, filename), dpi=300)
    #plt.show()

def plot_junction_scatter(df, library_name):
    """Standard Dot-Plot: Start vs Stop coordinates."""
    subset = df[df['Library'] == library_name]
    
    if subset.empty:
        return

    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Scatter plot with size proportional to count
    scatter = ax.scatter(subset['Start'], subset['Stop'], 
                         s=subset['Count'], alpha=0.6, c='purple', edgecolors='black')
    
    # Diagonal line (where Start == Stop)
    limit = max(subset[['Start', 'Stop']].max().max(), GENOME_LENGTH)
    ax.plot([0, limit], [0, limit], ls='--', c='gray', alpha=0.5)
    
    ax.set_xlim(0, limit)
    ax.set_ylim(0, limit)
    ax.set_xlabel("Donor Site (Start)")
    ax.set_ylabel("Acceptor Site (Stop)")
    ax.set_title(f"Recombination Junction Map\n{library_name}")
    ax.grid(True, linestyle=':', alpha=0.6)
    
    # Legend for sizes
    # Create dummy handles for legend
    sizes = [10, 50, 100, 300]
    legend_handles = [plt.scatter([], [], s=s, c='purple', edgecolors='black', alpha=0.6) for s in sizes]
    ax.legend(legend_handles, [f"{s} reads" for s in sizes], title="Frequency", loc='upper left')

    plt.tight_layout()
    filename = f"scatter_{_safe_filename(library_name)}.png"
    fig.savefig(os.path.join(OUTPUT_DIR, filename), dpi=300)
    #plt.show()

def plot_top_events(df, library_name):
    """Bar chart of the top 10 most frequent events."""
    subset = df[df['Library'] == library_name].copy()
    
    if subset.empty:
        return

    # Create a label column
    subset['Label'] = subset['Start'].astype(str) + " -> " + subset['Stop'].astype(str)
    
    # Get top 15
    top_events = subset.nlargest(15, 'Count')
    
    fig, ax = plt.subplots(figsize=(12, 6))
    bars = ax.bar(top_events['Label'], top_events['Count'], color='teal')
    
    ax.set_ylabel("Read Count")
    ax.set_xlabel("Junction Coordinates (Start -> Stop)")
    ax.set_title(f"Top 15 Recombination Events\n{library_name}")
    plt.xticks(rotation=45, ha='right')
    
    # Add count labels on top of bars
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height}', xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom')

    plt.tight_layout()
    filename = f"top_events_{_safe_filename(library_name)}.png"
    fig.savefig(os.path.join(OUTPUT_DIR, filename), dpi=300)
    #plt.show()

# --- Main Execution ---
if __name__ == "__main__":
    print(f"Reading file: {INPUT_FILE}...")
    df = parse_virema_results(INPUT_FILE)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    if not df.empty:
        print("Data loaded successfully.")
        print("Available Libraries in your file:")
        print(df['Library'].unique())
        
        # You can overwrite TARGET_LIBRARY here if you see a different name in the print output above
        print(f"\nVisualizing: {TARGET_LIBRARY}")
        
        plot_arc_diagram(df, TARGET_LIBRARY)
        plot_junction_scatter(df, TARGET_LIBRARY)
        plot_top_events(df, TARGET_LIBRARY)
    else:
        print("Error: DataFrame is empty. Check input file path and format.")
