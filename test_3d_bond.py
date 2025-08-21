# app_bond_angles_labeled_atoms.py
import pickle
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import streamlit as st

# --- Helper Functions ---
def save_molecule(mol, filename="molecule.pkl"):
    with open(filename, "wb") as file:
        pickle.dump(mol, file)

def load_molecule(filename="molecule.pkl"):
    with open(filename, "rb") as file:
        mol = pickle.load(file)
    return mol

def calculate_bond_angles(mol):
    bond_angles = []
    conf = mol.GetConformer()
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        neighbors1 = [a.GetIdx() for a in mol.GetAtomWithIdx(a1).GetNeighbors() if a.GetIdx() != a2]
        neighbors2 = [a.GetIdx() for a in mol.GetAtomWithIdx(a2).GetNeighbors() if a.GetIdx() != a1]
        for n1 in neighbors1:
            for n2 in neighbors2:
                pos1 = np.array(conf.GetAtomPosition(n1))
                pos2 = np.array(conf.GetAtomPosition(a1))
                pos3 = np.array(conf.GetAtomPosition(a2))
                vec1 = pos1 - pos2
                vec2 = pos3 - pos2
                angle = np.degrees(np.arccos(np.dot(vec1, vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))))
                bond_angles.append((a1,a2,angle))
    return bond_angles

def visualize_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        st.error("❌ Invalid SMILES string.")
        return None, None

    # Generate 3D coordinates
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol)
    save_molecule(mol)

    bond_angles = calculate_bond_angles(mol)
    conf = mol.GetConformer()

    view = py3Dmol.view(width=700, height=500)
    block = Chem.MolToMolBlock(mol)
    view.addModel(block, "mol")

    # --- Global stick style for crisp bonds ---
    view.setStyle({"stick": {"radius":0.2}})

    # --- Labeled spheres for atoms ---
    atom_colors = {"C":"gray","H":"white","O":"red","N":"blue",
                   "Cl":"green","Br":"brown","F":"cyan","P":"orange","S":"yellow"}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        symbol = atom.GetSymbol()
        color = atom_colors.get(symbol,"pink")
        pos = conf.GetAtomPosition(idx)
        
        # Draw sphere
        view.addSphere({"center":{"x":pos.x,"y":pos.y,"z":pos.z},"radius":0.35,"color":color})
        
        # Add label on sphere
        label_text = f"{symbol}({idx})"
        view.addLabel(label_text,
                      {"position":{"x":pos.x,"y":pos.y,"z":pos.z},
                       "backgroundColor":"black",
                       "fontColor":"white",
                       "fontSize":12,
                       "borderThickness":1,
                       "borderColor":"white",
                       "padding":2})

    # --- Bond angles as labels at bond midpoints ---
    for a1,a2,angle in bond_angles:
        pos1 = conf.GetAtomPosition(a1)
        pos2 = conf.GetAtomPosition(a2)
        mid_x = (pos1.x + pos2.x)/2
        mid_y = (pos1.y + pos2.y)/2
        mid_z = (pos1.z + pos2.z)/2
        view.addLabel(f"{angle:.1f}°",
                      {"position":{"x":mid_x,"y":mid_y,"z":mid_z},
                       "backgroundColor":"lightblue",
                       "fontColor":"black",
                       "fontSize":12,
                       "borderThickness":1})

    view.setBackgroundColor("white")
    view.zoomTo()
    view.spin(True)

    return view._make_html(), bond_angles

# --- Streamlit UI ---
st.set_page_config(page_title="Molecule Bond Angle Viewer", layout="wide")
st.title("🔬 Molecule Bond Angle Viewer")

smiles_input = st.text_input("Enter SMILES string", value="C1=CC=CC=C1")

if smiles_input:
    html, bond_angles = visualize_molecule(smiles_input)
    if html:
        st.components.v1.html(html, height=500, width=700, scrolling=False)
        st.subheader("📐 Bond Angles")
        st.dataframe([{"Atom1":a1,"Atom2":a2,"Angle (°)":round(angle,2)} for a1,a2,angle in bond_angles])

if st.button("Load Last Saved Molecule"):
    try:
        mol = load_molecule()
        block = Chem.MolToMolBlock(mol)
        view = py3Dmol.view(width=700, height=500)
        view.addModel(block, "mol")
        view.setStyle({"stick":{"radius":0.2}})
        view.setBackgroundColor("white")
        view.zoomTo()
        st.components.v1.html(view._make_html(), height=500, width=700, scrolling=False)
    except FileNotFoundError:
        st.error("⚠️ No saved molecule found yet.")
