"""
Module for 3D molecular visualization using py3Dmol and RDKit.
"""

import numpy as np
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
from src.utils.logging_utils import setup_logger

logger = setup_logger(__name__)

class Molecule3DVisualizer:
    """
    Class for creating 3D molecular visualizations.
    """
    
    def __init__(self):
        """Initialize the 3D visualizer."""
        logger.info("3D Molecule visualizer initialized")
    
    def save_molecule(self, mol, filename="molecule.pkl"):
        """Save RDKit molecule as a Pickle file."""
        try:
            with open(filename, "wb") as file:
                pickle.dump(mol, file)
            logger.info(f"Molecule saved as {filename}")
        except Exception as e:
            logger.error(f"Error saving molecule: {str(e)}")
    
    def load_molecule(self, filename="molecule.pkl"):
        """Load RDKit molecule from a Pickle file."""
        try:
            with open(filename, "rb") as file:
                mol = pickle.load(file)
            logger.info(f"Molecule loaded from {filename}")
            return mol
        except Exception as e:
            logger.error(f"Error loading molecule: {str(e)}")
            return None
    
    def calculate_bond_angles(self, mol):
        """Calculate bond angles of the molecule."""
        try:
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
                        bond_angles.append((n1, a1, a2, n2, angle))
            
            logger.info(f"Calculated {len(bond_angles)} bond angles")
            return bond_angles
        except Exception as e:
            logger.error(f"Error calculating bond angles: {str(e)}")
            return []
    
    def visualize_molecule_standard(self, smiles, chemical_name="Molecule"):
        """Generate standard 3D visualization - returns HTML for Streamlit."""
        try:
            mol = Chem.MolFromSmiles(smiles)
            
            if mol is None:
                logger.error("Invalid SMILES string")
                return None, []

            # Generate 3D coordinates
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            AllChem.MMFFOptimizeMolecule(mol)

            # Save the molecule
            self.save_molecule(mol, "saved_molecule.pkl")

            # Calculate bond angles
            bond_angles = self.calculate_bond_angles(mol)

            # Create 3Dmol viewer
            view = py3Dmol.view(width=800, height=600)
            block = Chem.MolToMolBlock(mol)
            view.addModel(block, "mol")
            view.setStyle({"stick": {}})
            view.setBackgroundColor("white")
            view.zoomTo()
            
            logger.info(f"Generated standard 3D visualization for {chemical_name}")
            return view, bond_angles
        
        except Exception as e:
            logger.error(f"Error in standard visualization: {str(e)}")
            return None, []
    
    def visualize_molecule_enhanced(self, smiles, chemical_name="Molecule"):
        """Generate enhanced 3D visualization with labeled atoms and bond angles."""
        try:
            mol = Chem.MolFromSmiles(smiles)
            
            if mol is None:
                logger.error("Invalid SMILES string")
                return None

            # Generate 3D coordinates
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            AllChem.MMFFOptimizeMolecule(mol)

            # Save the molecule
            self.save_molecule(mol, "saved_molecule.pkl")

            # Calculate bond angles
            bond_angles = self.calculate_bond_angles(mol)
            conf = mol.GetConformer()

            # Create 3Dmol viewer
            view = py3Dmol.view(width=800, height=600)
            block = Chem.MolToMolBlock(mol)
            view.addModel(block, "mol")

            # Global stick style for crisp bonds
            view.setStyle({"stick": {"radius": 0.2}})

            # Labeled spheres for atoms
            atom_colors = {
                "C": "gray", "H": "white", "O": "red", "N": "blue",
                "Cl": "green", "Br": "brown", "F": "cyan", "P": "orange", "S": "yellow"
            }
            
            for atom in mol.GetAtoms():
                idx = atom.GetIdx()
                symbol = atom.GetSymbol()
                color = atom_colors.get(symbol, "pink")
                pos = conf.GetAtomPosition(idx)
                
                # Draw sphere
                view.addSphere({
                    "center": {"x": pos.x, "y": pos.y, "z": pos.z},
                    "radius": 0.35,
                    "color": color
                })
                
                # Add label on sphere
                label_text = f"{symbol}({idx})"
                view.addLabel(label_text, {
                    "position": {"x": pos.x, "y": pos.y, "z": pos.z},
                    "backgroundColor": "black",
                    "fontColor": "white",
                    "fontSize": 12,
                    "borderThickness": 1,
                    "borderColor": "white",
                    "padding": 2
                })

            # Bond angles as labels at bond midpoints (limit to first 10 to avoid clutter)
            for i, (n1, a1, a2, n2, angle) in enumerate(bond_angles[:10]):
                pos1 = conf.GetAtomPosition(a1)
                pos2 = conf.GetAtomPosition(a2)
                mid_x = (pos1.x + pos2.x) / 2
                mid_y = (pos1.y + pos2.y) / 2
                mid_z = (pos1.z + pos2.z) / 2
                view.addLabel(f"{angle:.1f}°", {
                    "position": {"x": mid_x, "y": mid_y, "z": mid_z},
                    "backgroundColor": "lightblue",
                    "fontColor": "black",
                    "fontSize": 12,
                    "borderThickness": 1
                })

            view.setBackgroundColor("white")
            view.zoomTo()
            view.spin(True)

            logger.info(f"Generated enhanced 3D visualization with labels for {chemical_name}")
            return view
        
        except Exception as e:
            logger.error(f"Error in enhanced visualization: {str(e)}")
            return None
