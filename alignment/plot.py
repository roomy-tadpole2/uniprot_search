import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from Bio import AlignIO

def hex_to_rgb(hex_color):
    """Convert hex color (e.g., '#f94144') to an RGB tuple (0-1 range)."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) / 255 for i in (0, 2, 4))

def get_best_text_color(bg_color_hex):
    """Return 'black' or 'white' based on background brightness."""
    # Convert hex to RGB
    r, g, b = hex_to_rgb(bg_color_hex)
    
    # Calculate luminance
    luminance = 0.299 * r + 0.587 * g + 0.114 * b
    
    # Choose text color based on luminance
    return 'black' if luminance > 0.5 else 'white'

def color_adjust(color: str):
    return get_best_text_color(color)

def ylabel_color_map(name: str):
    if ('H2A' in name): return '#CC978E'
    if ('HTB' in name): return '#F39C6B'
    if ('HTR' in name): return '#FF3864'
    return 'black'

def plot(input_file_path, output_file_path: str, name_title=False, y_zip=1, ylabel_color=False):
    # Define colors for amino acids based on ESPript-like scheme
    color_map1 = {
        "Hydrophobic": "#f94144",
        "Aromatic": "#f8961e",
        "Negative": "#90be6d",
        "Positive": "#4d908e",
        "Polar uncharged": "#e15af1",
        "Sulfur_containing": "#ffe640",
        "Special": "#a186be"
    }
    color_map2 = {
        "Hydrophobic": "#ffd7e6",
        "Aromatic": "#ffd7e6",
        "Negative": "#d3f4ec",
        "Positive": "#d3f4ec",
        "Polar uncharged": "#fbded3",
        "Sulfur_containing": "#fbf3d9",
        "Special": "#a186be"
    }
    color_map = color_map1
    aa_color_scheme = {
        "A": color_map['Hydrophobic'], 
        "V": color_map['Hydrophobic'], 
        "L": color_map['Hydrophobic'], 
        "I": color_map['Hydrophobic'],  # Hydrophobic
        "F": color_map['Aromatic'], 
        "Y": color_map['Aromatic'], 
        "W": color_map['Aromatic'],    # Aromatic
        "D": color_map['Negative'], 
        "E": color_map['Negative'],    # Negative (acidic)
        "R": color_map['Positive'], 
        "K": color_map['Positive'], 
        "H": color_map['Positive'],    # Positive (basic)
        "S": color_map['Polar uncharged'], 
        "T": color_map['Polar uncharged'], 
        "N": color_map['Polar uncharged'], 
        "Q": color_map['Polar uncharged'], # Polar uncharged
        "C": color_map['Sulfur_containing'], 
        "M": color_map['Sulfur_containing'], # Sulfur-containing
        "G": '#a3a3a3', 
        "P": color_map['Special'],     # Special cases
        "-": "white"                   # Gaps
    }
    # Load aligned sequences
    alignment = AlignIO.read(f"{input_file_path}", "fasta")

    # Convert sequences to a matrix (rows = sequences, columns = positions)
    seq_data = [list(record.seq) for record in alignment]
    seq_labels = [record.id for record in alignment]
    num_seqs = len(seq_data)
    num_positions = len(seq_data[0])

    fig, ax = plt.subplots(figsize=(num_positions * 0.2, num_seqs * 0.5 * y_zip))

    for i, seq in enumerate(seq_data):
        for j, aa in enumerate(seq):
            color = aa_color_scheme.get(aa, "#0a0a0a")  # Default to black if unknown
            if (aa!='-'):
                ax.text(j, i*y_zip, aa, fontdict={'fontsize': 10, 'weight': 'heavy'}, 
                        color=color_adjust(color),
                        ha="center", va="center", bbox=
                        dict(facecolor=color, edgecolor=color, boxstyle='round,pad=0.2'))
            else:
                ax.text(j, i*y_zip, aa, fontsize=20, color="black",
                        ha="center", va="center", bbox=None)

    ax.set_xticks(range(num_positions))
    ax.set_xticklabels(range(1, num_positions + 1), fontsize=8, rotation=90)
    ax.set_yticks(np.array(range(num_seqs))*y_zip)

    ax.set_yticklabels(seq_labels, fontdict={'fontsize': 14, 'weight': 'heavy'})

    if (ylabel_color==True):
        print("setting y_label_colors...")
        for tick_label in ax.get_yticklabels():
            tick_label.set_color(ylabel_color_map(tick_label.get_text()))
    else:
        print("using defual black color for y_labels")
    
    ax.set_xlim(-0.5, num_positions - 0.5)

    print(f"num seq = {num_seqs}")

    ax.set_ylim(num_seqs*y_zip - 0.5*y_zip, -0.5*y_zip)
    ax.grid(False)

    # ax.yaxis.set_major_locator(ticker.MultipleLocator(0.8))
    yticks = ax.get_yticks()
    spacing = yticks[1] - yticks[0]  # Compute spacing
    print(f"Current y-tick spacing: {spacing}")

    # Remove axes
    ax.set_frame_on(False)
    ax.xaxis.set_tick_params(size=0)
    ax.yaxis.set_tick_params(size=0)

    plt.tight_layout()

    if (name_title):
        plt.title(output_file_path, fontsize=35, fontfamily='serif', color='#8B4513', loc='left', pad=20)  

    if (output_file_path!=None):
        plt.savefig(f"{output_file_path}", dpi=300, bbox_inches="tight")
    else:plt.show()
