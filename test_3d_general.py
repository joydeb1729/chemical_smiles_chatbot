# app.py
import pickle
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
import streamlit as st

# --- Helper Functions ---
def save_molecule(mol, filename="molecule.pkl"):
    with open(filename, "wb") as file:
        pickle.dump(mol, file)

def load_molecule(filename="molecule.pkl"):
    with open(filename, "rb") as file:
        mol = pickle.load(file)
    return mol

def visualize_molecule(smiles):
    """Return 3Dmol HTML for a molecule from SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol)
    save_molecule(mol)

    block = Chem.MolToMolBlock(mol)
    view = py3Dmol.view(width=500, height=400)
    view.addModel(block, "mol")
    view.setStyle({"stick": {}})
    view.setBackgroundColor("white")
    view.zoomTo()
    return view._make_html()

# --- Streamlit UI ---
st.set_page_config(page_title="3D Molecule Viewer", layout="centered")
st.title("🔬 3D Molecule Visualizer")

# ✅ Input box always visible
smiles_input = st.text_input(
    "Enter a SMILES string", 
    value="C1=CC=CC=C1",  # default benzene
    placeholder="e.g. C(Cl)Cl for dichloromethane"
)

if smiles_input:
    html = visualize_molecule(smiles_input)
    if html:
        st.components.v1.html(html, height=500, width=700, scrolling=False)
    else:
        st.error("Invalid SMILES string.")

if st.button("Load Last Saved Molecule"):
    try:
        mol = load_molecule()
        block = Chem.MolToMolBlock(mol)
        view = py3Dmol.view(width=500, height=400)
        view.addModel(block, "mol")
        view.setStyle({"stick": {}})
        view.zoomTo()
        st.components.v1.html(view._make_html(), height=500, width=700, scrolling=False)
    except FileNotFoundError:
        st.error("⚠️ No saved molecule found yet.")
