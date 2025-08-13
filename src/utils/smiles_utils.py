"""
Utility functions for SMILES and chemical data processing.
"""

import re

def clean_chemical_name(name):
    """
    Clean and normalize a chemical name.
    
    Args:
        name (str): The chemical name to clean.
        
    Returns:
        str: The cleaned chemical name.
    """
    if not name:
        return name
    
    # Remove any leading/trailing whitespace
    name = name.strip()
    
    # Remove any parentheses and their contents
    name = re.sub(r'\s*\([^)]*\)', '', name)
    
    # Capitalize the first letter of each word
    name = ' '.join(word.capitalize() for word in name.split())
    
    return name

def validate_smiles(smiles):
    """
    Validate a SMILES string.
    
    This is a basic validation that checks for common issues.
    For a more thorough validation, use RDKit's parsing.
    
    Args:
        smiles (str): The SMILES string to validate.
        
    Returns:
        bool: True if the SMILES string appears valid, False otherwise.
    """
    if not smiles:
        return False
    
    # Check for balanced parentheses and brackets
    if smiles.count('(') != smiles.count(')'):
        return False
    if smiles.count('[') != smiles.count(']'):
        return False
    
    # Check for common atoms
    common_atoms = ['C', 'H', 'O', 'N', 'P', 'S', 'F', 'Cl', 'Br', 'I']
    has_atom = False
    
    for atom in common_atoms:
        if atom in smiles:
            has_atom = True
            break
    
    if not has_atom:
        # If no common atoms found, check for other valid elements
        element_pattern = r'\[([A-Z][a-z]?)(?:[^]]*)]'
        elements = re.findall(element_pattern, smiles)
        if not elements:
            return False
    
    return True
