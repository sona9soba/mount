import pandas as pd
import pymol
import sys
import os

COLORS = ["red", "blue", "green", "yellow", "magenta", "cyan", "orange", "purple", 
          "lime", "pink", "teal", "lavender", "brown", "beige", "maroon", "mint", 
          "olive", "apricot", "navy", "grey"] 
          

          
def generate_dicts_from_csv(file_path):
    df = pd.read_csv(file_path)
    grouped = df.groupby('uniprot_id')
    
    uniprot_ids = []
    dicts_list = []
    for uniprot_id, group in grouped:
        uniprot_ids.append(uniprot_id)
        rank_residue_mapping = {}
        for _, row in group.iterrows():
            rank_residue_mapping[row['rank']] = [residue.strip() for residue in row['residue_id'].split(',')]
        dicts_list.append(rank_residue_mapping)
    return uniprot_ids, dicts_list
    


def color_by_group(pdb_filename, residue_dict):
    
    pymol = pymol.PyMol()
    pymol.finish_launching(['pymol', '-cQi'])  # Launch PyMOL in command-line and quiet mode

    # Step 1: Load the PDB file
    pymol.cmd.load(pdb_filename, "protein")

    # Step 2: Show in surface mode
    pymol.cmd.show("surface", "protein")

    # Step 3 & 4: Color code by groups
    for key in sorted(residue_dict.keys()):
        color = COLORS[key % len(COLORS)]  # Cycle through colors if there are more groups than colors
        # Extract only the index numbers from the residue identifiers
        indices = [residue.lstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ") for residue in residue_dict[key]]
        residues = "+".join(indices)
        selection_string = f"protein and resi {residues}"
        pymol.cmd.color(color, selection_string)

    # Step 6: Set the view to focus on the highest-ranked residues
    top_group_indices = [residue.lstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ") for residue in residue_dict[min(residue_groups.keys())]]
    pymol.cmd.zoom(f"protein and resi {'+'.join(top_group_indices)}", buffer=2.0)

    # Export the file
    output_base = pdb_filename.replace(".pdb", "")
    pymol.cmd.png(f"{output_base}.png", dpi=300, height=2400)
    pymol.cmd.save(f"{output_base}.pse")

    # Quit PyMOL
    pymol.cmd.quit()


input_csv = 'modified.csv'
uniprot_list, dicts_list = generate_dicts_from_csv(input_csv)
pdb_list = [uniprot_id + '.pdb' for uniprot_id in uniprot_list]

print(len(uniprot_list), len(dicts_list), len(pdb_list))

for i in range(len(uniprot_list)):
    pdb_filename = pdb_list[i]
    residue_dict = dicts_list[i]
    color_by_group(pdb_filename, residue_dict)
   

  
