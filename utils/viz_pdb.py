# visualization utils for BeadsNet
#!pip install biopython py3Dmol

from Bio.PDB import PDBList, PDBParser, PDBIO, Select
import py3Dmol
import sys
from IPython.display import display
import ipywidgets as widgets
from ipywidgets import interact

# --- 1. Fetch the PDB Structure and Create C-alpha file ---
# This initial setup is required for the different visualizations.
pdb_id = '9NQU'
pdbl = PDBList()

try:
    # Fetch the full PDB file
    original_pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')
    print(f"PDB file for {pdb_id} saved to: {original_pdb_file}")

    # --- Create a C-alpha only PDB file ---
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, original_pdb_file)

    class CalphaSelect(Select):
        """A selector class to extract only C-alpha atoms."""
        def accept_chain(self, chain):
            return True
        def accept_residue(self, residue):
            return True
        def accept_atom(self, atom):
            return atom.get_name() == 'CA'

    calpha_pdb_file = f"{pdb_id}_calpha.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(calpha_pdb_file, CalphaSelect())
    #print(f"C-alpha atoms file created at: {calpha_pdb_file}")

except Exception as e:
    print(f"An error occurred during file setup: {e}")
    sys.exit()


# --- 2. Read File Contents for Visualization ---
# We read the file contents into strings to pass them to py3Dmol.
try:
    with open(original_pdb_file, "r") as f1:
        full_pdb_data = f1.read()
    with open(calpha_pdb_file, "r") as f2:
        calpha_pdb_data = f2.read()
except FileNotFoundError as e:
    print(f"Error reading PDB files: {e}")
    sys.exit()


# --- 3. Create Interactive Visualization ---
#print("\n## Interactive Protein Visualizer")

# Define the function that creates the view based on the selected style.
# This function will be called by 'interact' every time the dropdown changes.
def create_interactive_view(style):
    # Create a new viewer instance each time
    view = py3Dmol.view(width=800, height=800)

    # Add both models (full protein and C-alpha only)
    view.addModel(full_pdb_data, "pdb")
    view.addModel(calpha_pdb_data, "pdb")

    # Apply the selected style
    if style == 'Cartoon':
        # Show only the full model (model 0)
        view.setStyle({'model': 0}, {'cartoon': {'color': 'spectrum'}})
        view.setStyle({'model': 1}, {}) # Hide C-alpha model

    elif style == 'Cartoon with Surface':
        # 1) Draw the cartoon + sticks as before
        view.setStyle({'model': 0}, {'cartoon': {'color': 'spectrum'}})
        view.setStyle({'model': 1}, {}) # Hide C-alpha model
        # 2) Surface only the 20 amino‑acid residues in model 0
        protein_residues = [
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
            "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
        ]
        view.addSurface(
            py3Dmol.VDW,
            {'opacity': 0.8, 'color': 'white'},
            {'model': 0, 'resn': protein_residues}
        )
        # 3) Hide your Cα‐only model
        view.setStyle({'model': 1}, {})

    elif style == 'Cartoon with Residue':
        view.setStyle({'model': 0}, {'cartoon': {'color': 'spectrum'}})
        view.addStyle({'model': 0}, {'stick': {}})
        view.setStyle({'model': 1}, {})
    elif style == 'Cartoon with Residue + Surface':
        # 1) Draw the cartoon + sticks as before
        view.setStyle({'model': 0}, {'cartoon': {'color': 'spectrum'}})
        view.addStyle({'model': 0}, {'stick': {}})
        # 2) Surface only the 20 amino‑acid residues in model 0
        protein_residues = [
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
            "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
        ]
        view.addSurface(
            py3Dmol.VDW,
            {'opacity': 0.8, 'color': 'white'},
            {'model': 0, 'resn': protein_residues}
        )
        # 3) Hide your Cα‐only model
        view.setStyle({'model': 1}, {})

    elif style == 'Stick and Ball':
        view.setStyle({'model': 0}, {'stick': {}})
        view.setStyle({'model': 1}, {})
    elif style == 'C-alpha Backbone Trace':
        # Show only the C-alpha model (model 1)
        view.setStyle({'model': 1}, {'cartoon': {'style': 'trace', 'color': 'spectrum', 'thickness': 0.1}})
        view.setStyle({'model': 0}, {}) # Hide full model

    view.zoomTo()
    # Add the spin/turntable effect
    view.spin()
    # Return the view object to be displayed by interact
    return view.show()

# Create the dropdown widget for style selection
style_options = [
    'Cartoon',
    'Cartoon with Surface',    
    'Cartoon with Residue',
    'Cartoon with Residue + Surface',
    'Stick and Ball',
    'C-alpha Backbone Trace'
]

# Use the interact function to link the dropdown to our view creation function
interact(create_interactive_view, style=widgets.Dropdown(options=style_options, description='Select Style:'));
