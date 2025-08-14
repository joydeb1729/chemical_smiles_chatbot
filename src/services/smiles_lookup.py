"""
Module for looking up SMILES representations from chemical names.
"""

import pubchempy as pcp
import os
import json
import time
from datetime import datetime
from py2opsin import py2opsin

from src.utils.logging_utils import setup_logger
from config.settings import PUBCHEM_TIMEOUT, CACHE_DIR, CACHE_EXPIRY

logger = setup_logger(__name__)

class SMILESLookup:
    """
    Class to look up SMILES representations of chemicals using OPSIN and PubChem.
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
        
        # Debug: Log the exact chemical name and format we're looking up
        logger.info(f"Looking up chemical: '{chemical_name}' (normalized: '{chemical_name_lower}')")
        
        # Check cache first if enabled
        if self.use_cache and chemical_name_lower in self.cache:
            cache_entry = self.cache[chemical_name_lower]
            if self._is_cache_valid(cache_entry):
                logger.info(f"Found cached SMILES for {chemical_name}")
                return cache_entry['smiles']
            else:
                logger.info(f"Cached SMILES for {chemical_name} expired")
        
        # First try using py2opsin (a Python wrapper for OPSIN)
        try:
            logger.info(f"Looking up SMILES for {chemical_name} using py2opsin")
            
            # Try to get SMILES using py2opsin
            smiles = py2opsin(chemical_name=chemical_name, output_format="SMILES")
            
            # Check if we got a valid result
            if smiles and len(smiles) > 0:
                logger.info(f"Found SMILES for {chemical_name} using py2opsin: {smiles}")
                
                # Cache the result if caching is enabled
                if self.use_cache:
                    self.cache[chemical_name_lower] = {
                        'smiles': smiles,
                        'timestamp': time.time(),
                        'date': datetime.now().isoformat(),
                        'source': 'py2opsin'
                    }
                    self._save_cache()
                
                return smiles
            else:
                logger.warning(f"py2opsin returned empty result for {chemical_name}")
        except Exception as e:
            logger.warning(f"py2opsin lookup failed for {chemical_name}: {str(e)}")
        
        # Try with alternate chemical name formats
        alt_names = [
            chemical_name.replace(" ", ""),  # No spaces
            chemical_name.replace("-", ""),  # No hyphens
            chemical_name.replace(" ", "-")  # Replace spaces with hyphens
        ]
        
        for alt_name in alt_names:
            try:
                logger.info(f"Trying alternate format with py2opsin: {alt_name}")
                
                smiles = py2opsin(chemical_name=alt_name, output_format="SMILES")
                
                if smiles and len(smiles) > 0:
                    logger.info(f"Found SMILES for alternate name {alt_name} using py2opsin: {smiles}")
                    
                    # Cache the result if caching is enabled
                    if self.use_cache:
                        self.cache[chemical_name_lower] = {
                            'smiles': smiles,
                            'timestamp': time.time(),
                            'date': datetime.now().isoformat(),
                            'source': 'py2opsin_alternate'
                        }
                        self._save_cache()
                    
                    return smiles
            except Exception as e:
                logger.warning(f"py2opsin alternate format lookup failed for {alt_name}: {str(e)}")
        
        # If py2opsin fails, try PubChem
        try:
            logger.info(f"Looking up SMILES for {chemical_name} using PubChem")
            
            # Search for the compound in PubChem
            compounds = pcp.get_compounds(chemical_name, 'name', timeout=PUBCHEM_TIMEOUT)
            
            if compounds:
                # Get the SMILES for the first (most relevant) compound
                compound = compounds[0]
                
                # Verify we have a valid compound with a canonical SMILES
                if hasattr(compound, 'canonical_smiles') and compound.canonical_smiles:
                    smiles = compound.canonical_smiles
                    
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
                    logger.warning(f"Compound found for {chemical_name} in PubChem, but no SMILES available")
            else:
                logger.warning(f"No compounds found for {chemical_name} in PubChem")
        except Exception as e:
            logger.error(f"PubChem lookup failed for {chemical_name}: {str(e)}")
        
        logger.error(f"Could not find SMILES for {chemical_name} using any method")
        return None
