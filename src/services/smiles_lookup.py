"""
Module for looking up SMILES representations from chemical names.
"""

import pubchempy as pcp
import os
import json
import time
import requests
from datetime import datetime

from src.utils.logging_utils import setup_logger
from config.settings import PUBCHEM_TIMEOUT, CACHE_DIR, CACHE_EXPIRY

logger = setup_logger(__name__)

# Mapping of common chemicals with their SMILES representations
COMMON_CHEMICALS = {
    "benzene": "C1=CC=CC=C1",
    "water": "O",
    "methane": "C",
    "ethanol": "CCO",
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "acetaminophen": "CC(=O)NC1=CC=C(C=C1)O",
    "glucose": "C(C1C(C(C(C(O1)O)O)O)O)O",
    "sodium chloride": "[Na+].[Cl-]",
    "methylbenzene": "CC1=CC=CC=C1",  # toluene
    "ethylbenzene": "CCC1=CC=CC=C1",
    "propylbenzene": "CCCC1=CC=CC=C1",
    "carbon dioxide": "O=C=O",
    "benzaldehyde": "C1=CC=C(C=C1)C=O",
    "acetone": "CC(=O)C",
    "ethylene": "C=C",
    "propylene": "CC=C",
    "ammonia": "N",
    "hydrogen peroxide": "OO",
    "formaldehyde": "C=O",
    "acetic acid": "CC(=O)O"
}

class SMILESLookup:
    """
    Class to look up SMILES representations of chemicals using PubChem.
    """
    
    def __init__(self, use_cache=True):
        """
        Initialize the SMILES lookup service.
        
        Args:
            use_cache (bool, optional): Whether to use the cache. Defaults to True.
        """
        self.use_cache = use_cache
        self.cache_file = os.path.join(CACHE_DIR, "smiles_cache.json")
        self.cache = {}
        
        # Create cache directory if it doesn't exist
        if self.use_cache:
            os.makedirs(CACHE_DIR, exist_ok=True)
            self._load_cache()
    
    def _load_cache(self):
        """
        Load the cache from disk.
        """
        if os.path.exists(self.cache_file):
            try:
                with open(self.cache_file, 'r') as f:
                    self.cache = json.load(f)
                logger.info(f"Loaded SMILES cache with {len(self.cache)} entries")
            except Exception as e:
                logger.error(f"Error loading SMILES cache: {str(e)}")
                self.cache = {}
    
    def _save_cache(self):
        """
        Save the cache to disk.
        """
        try:
            with open(self.cache_file, 'w') as f:
                json.dump(self.cache, f, indent=2)
            logger.info(f"Saved SMILES cache with {len(self.cache)} entries")
        except Exception as e:
            logger.error(f"Error saving SMILES cache: {str(e)}")
    
    def _is_cache_valid(self, entry):
        """
        Check if a cache entry is still valid.
        
        Args:
            entry (dict): The cache entry.
            
        Returns:
            bool: True if valid, False if expired.
        """
        timestamp = entry.get('timestamp', 0)
        current_time = time.time()
        return (current_time - timestamp) < CACHE_EXPIRY
    
    def get_smiles(self, chemical_name):
        """
        Get the SMILES representation for a chemical name.
        
        Args:
            chemical_name (str): The name of the chemical.
            
        Returns:
            str: The SMILES representation, or None if not found.
        """
        if not chemical_name:
            logger.error("Empty chemical name provided")
            return None
            
        # Normalize the chemical name
        chemical_name_lower = chemical_name.lower()
        
        # Normalize some common variations
        chemical_variations = {
            "paracetamol": "acetaminophen",
            "salt": "sodium chloride",
            "table salt": "sodium chloride",
            "toluene": "methylbenzene",
            "carbon di oxide": "carbon dioxide"
        }
        
        if chemical_name_lower in chemical_variations:
            chemical_name_lower = chemical_variations[chemical_name_lower]
            logger.info(f"Normalized chemical name: {chemical_name} -> {chemical_name_lower}")
        
        # Direct lookup for common chemicals without API calls
        if chemical_name_lower in COMMON_CHEMICALS:
            logger.info(f"Direct lookup for {chemical_name} in predefined list")
            smiles = COMMON_CHEMICALS[chemical_name_lower]
            
            # Cache the result if caching is enabled
            if self.use_cache:
                self.cache[chemical_name_lower] = {
                    'smiles': smiles,
                    'timestamp': time.time(),
                    'date': datetime.now().isoformat(),
                    'source': 'predefined'
                }
                self._save_cache()
                
            return smiles
        
        # Check cache first if enabled
        if self.use_cache and chemical_name_lower in self.cache:
            cache_entry = self.cache[chemical_name_lower]
            if self._is_cache_valid(cache_entry):
                logger.info(f"Found cached SMILES for {chemical_name}")
                return cache_entry['smiles']
            else:
                logger.info(f"Cached SMILES for {chemical_name} expired")
        
        # First try OPSIN service (more reliable for systematic names)
        try:
            logger.info(f"Looking up SMILES for {chemical_name} using OPSIN")
            # URL encode the chemical name to handle spaces and special characters
            import urllib.parse
            encoded_name = urllib.parse.quote(chemical_name)
            opsin_url = f"https://opsin.ch.cam.ac.uk/opsin/{encoded_name}.smiles"
            
            response = requests.get(opsin_url, timeout=10)
            
            # OPSIN returns an error message when it can't parse a name
            if response.status_code == 200 and response.text and "FAILED" not in response.text and "ERROR" not in response.text:
                smiles = response.text.strip()
                logger.info(f"Found SMILES for {chemical_name} using OPSIN: {smiles}")
                
                # Cache the result if caching is enabled
                if self.use_cache:
                    self.cache[chemical_name_lower] = {
                        'smiles': smiles,
                        'timestamp': time.time(),
                        'date': datetime.now().isoformat(),
                        'source': 'opsin'
                    }
                    self._save_cache()
                
                return smiles
            else:
                logger.warning(f"OPSIN could not parse chemical name: {chemical_name}. Response: {response.text[:100]}")
        except Exception as e:
            logger.warning(f"OPSIN lookup failed for {chemical_name}: {str(e)}")
            
        # Try with alternate chemical name formats
        alt_names = [
            chemical_name,
            chemical_name.replace(" ", ""),  # No spaces
            chemical_name.replace("-", ""),  # No hyphens
            chemical_name.replace(" ", "-")  # Replace spaces with hyphens
        ]
        
        for alt_name in alt_names:
            if alt_name == chemical_name:
                continue  # Skip the original name as we already tried it
                
            try:
                logger.info(f"Trying alternate format for OPSIN: {alt_name}")
                import urllib.parse
                encoded_name = urllib.parse.quote(alt_name)
                opsin_url = f"https://opsin.ch.cam.ac.uk/opsin/{encoded_name}.smiles"
                
                response = requests.get(opsin_url, timeout=10)
                
                if response.status_code == 200 and response.text and "FAILED" not in response.text and "ERROR" not in response.text:
                    smiles = response.text.strip()
                    logger.info(f"Found SMILES for alternate name {alt_name} using OPSIN: {smiles}")
                    
                    # Cache the result if caching is enabled
                    if self.use_cache:
                        self.cache[chemical_name_lower] = {
                            'smiles': smiles,
                            'timestamp': time.time(),
                            'date': datetime.now().isoformat(),
                            'source': 'opsin_alternate'
                        }
                        self._save_cache()
                    
                    return smiles
            except Exception as e:
                logger.warning(f"OPSIN alternate format lookup failed for {alt_name}: {str(e)}")
        
        # If OPSIN fails, try PubChem
        try:
            logger.info(f"Looking up SMILES for {chemical_name} using PubChem")
            
            # Search for the compound in PubChem
            compounds = pcp.get_compounds(chemical_name, 'name', timeout=PUBCHEM_TIMEOUT)
            
            if compounds:
                # Get the SMILES for the first (most relevant) compound
                smiles = compounds[0].canonical_smiles
                
                # Cache the result if caching is enabled
                if self.use_cache:
                    self.cache[chemical_name_lower] = {
                        'smiles': smiles,
                        'timestamp': time.time(),
                        'date': datetime.now().isoformat(),
                        'source': 'pubchem'
                    }
                    self._save_cache()
                
                logger.info(f"Found SMILES for {chemical_name} using PubChem: {smiles}")
                return smiles
            else:
                logger.warning(f"No compounds found for {chemical_name} in PubChem")
        except Exception as e:
            logger.error(f"PubChem lookup failed for {chemical_name}: {str(e)}")
        
        # If all online services fail, check if we have a predefined SMILES
        if chemical_name_lower in COMMON_CHEMICALS:
            logger.info(f"Using predefined SMILES for {chemical_name}")
            smiles = COMMON_CHEMICALS[chemical_name_lower]
            
            # Cache the result if caching is enabled
            if self.use_cache:
                self.cache[chemical_name_lower] = {
                    'smiles': smiles,
                    'timestamp': time.time(),
                    'date': datetime.now().isoformat(),
                    'source': 'predefined'
                }
                self._save_cache()
                
            return smiles
            
        # Check if any similar named chemicals are available
        for common_name, smiles in COMMON_CHEMICALS.items():
            if common_name in chemical_name_lower or chemical_name_lower in common_name:
                logger.info(f"Using similar chemical SMILES for {chemical_name} -> {common_name}")
                return smiles
        
        logger.error(f"Could not find SMILES for {chemical_name} using any method")
        return None
