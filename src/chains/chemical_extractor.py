"""
Module for chemical name extraction.
"""

import re
from src.utils.logging_utils import setup_logger
from config.prompts import CHEMICAL_EXTRACTION_TEMPLATE, CHEMICAL_INFORMATION_TEMPLATE

logger = setup_logger(__name__)

# Common chemicals that can be directly identified from queries
COMMON_CHEMICALS = {
    "benzene": ["benzene"],
    "water": ["water", "h2o", "h₂o"],
    "sodium chloride": ["salt", "sodium chloride", "nacl", "table salt", "kitchen salt"],
    "caffeine": ["caffeine"],
    "aspirin": ["aspirin", "acetylsalicylic acid"],
    "ethanol": ["ethanol", "alcohol", "ethyl alcohol"],
    "glucose": ["glucose", "dextrose"],
    "methane": ["methane", "ch4", "ch₄", "methain", "methan"],  # Added misspellings
    "acetaminophen": ["acetaminophen", "paracetamol", "tylenol"],
    "carbon dioxide": ["carbon dioxide", "co2", "co₂", "carbon di oxide", "carbondioxide"],
    "methylbenzene": ["methylbenzene", "toluene"],
    "ethylbenzene": ["ethylbenzene"],
    "propylbenzene": ["propylbenzene"],
    "benzaldehyde": ["benzaldehyde"],
    "acetone": ["acetone", "dimethyl ketone", "propanone"],
    "ethylene": ["ethylene", "ethene"],
    "propylene": ["propylene", "propene"],
    "ammonia": ["ammonia", "nh3", "nh₃"],
    "hydrogen peroxide": ["hydrogen peroxide", "h2o2", "h₂o₂"],
    "formaldehyde": ["formaldehyde", "methanal"],
    "acetic acid": ["acetic acid", "ethanoic acid", "vinegar acid"]
}

class ChemicalExtractor:
    """
    Class to extract chemical names from natural language queries.
    """
    
    def __init__(self, llm):
        """
        Initialize the ChemicalExtractor.
        
        Args:
            llm: The language model to use for extraction.
        """
        self.llm = llm
        logger.info("Chemical extractor initialized")
    
    def extract_chemical_name(self, query):
        """
        Extract a chemical name from a natural language query.
        
        Args:
            query (str): The natural language query.
            
        Returns:
            str: The extracted chemical name, or None if no chemical is found.
        """
        try:
            logger.info(f"Extracting chemical name from query: {query}")
            
            # First try to match common chemicals directly
            query_lower = query.lower()
            
            # Check for direct matches or common patterns like "what is [chemical]?"
            for chemical, aliases in COMMON_CHEMICALS.items():
                for alias in aliases:
                    # Check for exact mention
                    if alias in query_lower:
                        logger.info(f"Direct match found for chemical: {chemical}")
                        return chemical
                    
                    # Check for "what is X" pattern
                    what_is_pattern = f"what is {alias}"
                    if what_is_pattern in query_lower:
                        logger.info(f"'What is' pattern match found for chemical: {chemical}")
                        return chemical
            
            # Use the LLM with direct prompt for more flexibility
            full_prompt = CHEMICAL_EXTRACTION_TEMPLATE.format(query=query)
            logger.info(f"Using LLM with prompt: {full_prompt[:100]}...")
            
            # Use the LLM directly instead of through LangChain
            response = self.llm.generate(full_prompt)
            
            # Clean up the result
            chemical = response.strip()
            if chemical.lower() == "none" or not chemical:
                return None
            
            # Normalize common chemical variations
            normalized_chemical = self._normalize_chemical_name(chemical)
            if normalized_chemical:
                logger.info(f"Normalized chemical name: {chemical} -> {normalized_chemical}")
                return normalized_chemical
            
            logger.info(f"Extracted chemical name: {chemical}")
            return chemical
        
        except Exception as e:
            logger.error(f"Error extracting chemical name: {str(e)}")
            return None
    
    def _normalize_chemical_name(self, chemical_name):
        """
        Normalize chemical name variations to their standard form.
        
        Args:
            chemical_name (str): The chemical name to normalize.
            
        Returns:
            str: The normalized chemical name, or the original if no normalization is needed.
        """
        if not chemical_name:
            return None
            
        chemical_lower = chemical_name.lower()
        
        # Common variations mapping
        variations = {
            # Spacing variations
            "carbon di oxide": "carbon dioxide",
            "carbondioxide": "carbon dioxide",
            "carbon-di-oxide": "carbon dioxide",
            "sodiumchloride": "sodium chloride",
            "sodium-chloride": "sodium chloride",
            "hydrogenperoxide": "hydrogen peroxide",
            "hydrogen-peroxide": "hydrogen peroxide",
            # Common abbreviations
            "h2o": "water",
            "co2": "carbon dioxide",
            "nacl": "sodium chloride",
            "nh3": "ammonia",
            "h2o2": "hydrogen peroxide",
            # Common name variations
            "table salt": "sodium chloride",
            "kitchen salt": "sodium chloride",
            "vinegar acid": "acetic acid",
            "methylbenzene": "toluene",
            "paracetamol": "acetaminophen",
            "acetylsalicylic acid": "aspirin",
            "ethyl alcohol": "ethanol",
            "drinking alcohol": "ethanol"
        }
        
        # Check if the chemical name has a known variation
        if chemical_lower in variations:
            return variations[chemical_lower]
        
        # Try to remove spaces and check again
        no_spaces = chemical_lower.replace(" ", "")
        if no_spaces in variations:
            return variations[no_spaces]
            
        # Try with different space/hyphen formats
        hyphenated = chemical_lower.replace(" ", "-")
        if hyphenated in variations:
            return variations[hyphenated]
            
        # No normalization needed
        return None
    
    def get_chemical_information(self, chemical, smiles):
        """
        Get information about a chemical compound.
        
        Args:
            chemical (str): The name of the chemical.
            smiles (str): The SMILES representation of the chemical.
            
        Returns:
            str: Information about the chemical.
        """
        try:
            logger.info(f"Getting information for chemical: {chemical}")
            
            # Use the LLM directly instead of through LangChain
            prompt = CHEMICAL_INFORMATION_TEMPLATE.format(chemical=chemical, smiles=smiles)
            result = self.llm.generate(prompt)
            
            logger.info(f"Generated information for {chemical}")
            return result.strip()
        
        except Exception as e:
            logger.error(f"Error getting chemical information: {str(e)}")
            return f"Unable to retrieve information for {chemical}: {str(e)}"
    
    def get_chemical_information(self, chemical, smiles):
        """
        Get information about a chemical compound.
        
        Args:
            chemical (str): The name of the chemical.
            smiles (str): The SMILES representation of the chemical.
            
        Returns:
            str: Information about the chemical.
        """
        try:
            logger.info(f"Getting information for chemical: {chemical}")
            
            # Run the information chain
            result = self.info_chain.run(chemical=chemical, smiles=smiles)
            
            logger.info(f"Generated information for {chemical}")
            return result.strip()
        
        except Exception as e:
            logger.error(f"Error getting chemical information: {str(e)}")
            return f"Unable to retrieve information for {chemical}: {str(e)}"
