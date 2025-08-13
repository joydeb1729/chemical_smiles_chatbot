"""
Module for chemical name extraction.
"""

import re
from src.utils.logging_utils import setup_logger
from config.prompts import CHEMICAL_EXTRACTION_TEMPLATE, CHEMICAL_INFORMATION_TEMPLATE

logger = setup_logger(__name__)

class ChemicalExtractor:
    """
    Class to extract chemical names from natural language queries using LLM.
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
        Extract a chemical name from a natural language query using LLM.
        
        Args:
            query (str): The natural language query.
            
        Returns:
            str: The extracted chemical name, or None if no chemical is found.
        """
        try:
            logger.info(f"Extracting chemical name from query: {query}")
            
            # Use the LLM with the chemical extraction prompt
            full_prompt = CHEMICAL_EXTRACTION_TEMPLATE.format(query=query)
            logger.info(f"Using LLM with prompt: {full_prompt[:100]}...")
            
            # Use the LLM directly
            response = self.llm.generate(full_prompt)
            logger.debug(f"Raw LLM response: {response[:100]}")
            
            # Clean up the result
            chemical = response.strip()
            
            # Handle model returning the instruction text instead of just the chemical name
            if "provide the standardized chemical name" in chemical.lower():
                logger.warning("Model returned instruction text instead of chemical name")
                # Extract the chemical from the query directly
                query_lower = query.lower()
                
                # Try to extract using regex patterns for common query forms
                if "what is" in query_lower:
                    # Extract text after "what is"
                    chemical = query_lower.replace("what is", "").strip()
                    logger.info(f"Extracted from 'what is' pattern: {chemical}")
                elif "tell me about" in query_lower:
                    # Extract text after "tell me about"
                    chemical = query_lower.replace("tell me about", "").strip()
                    logger.info(f"Extracted from 'tell me about' pattern: {chemical}")
                elif "explain" in query_lower:
                    # Extract text after "explain"
                    chemical = query_lower.replace("explain", "").strip()
                    logger.info(f"Extracted from 'explain' pattern: {chemical}")
                elif re.search(r'of\s+([a-z0-9]+)', query_lower):
                    # Try to extract patterns like "properties of X"
                    match = re.search(r'of\s+([a-z0-9]+)', query_lower)
                    chemical = match.group(1).strip()
                    logger.info(f"Extracted from 'of X' pattern: {chemical}")
            
            # Final validation - if chemical is too long or contains instructional text
            if chemical and (len(chemical) > 30 or "chemical name" in chemical.lower() or 
                            "standardized" in chemical.lower() or "provide" in chemical.lower()):
                logger.warning(f"Invalid chemical name detected: {chemical}")
                # Extract from query as last resort
                query_words = query.lower().split()
                if len(query_words) > 0:
                    # Typically the last word in a "what is X" query is the chemical
                    chemical = query_words[-1].strip("?.,;:")
                    logger.info(f"Using last word from query as chemical: {chemical}")
            
            if not chemical or chemical.lower() == "none":
                logger.warning(f"No chemical detected in query: {query}")
                return None
            
            logger.info(f"Final extracted chemical name: {chemical}")
            return chemical
        
        except Exception as e:
            logger.error(f"Error extracting chemical name: {str(e)}")
            # Emergency fallback - try simple pattern matching on the query
            try:
                # Last resort extraction directly from query
                query_lower = query.lower()
                if "what is" in query_lower:
                    chemical = query_lower.replace("what is", "").strip()
                    logger.info(f"Emergency fallback extraction: {chemical}")
                    return chemical
            except:
                pass
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
            
            # Use the LLM directly with the chemical information prompt
            prompt = CHEMICAL_INFORMATION_TEMPLATE.format(chemical=chemical, smiles=smiles)
            result = self.llm.generate(prompt)
            
            logger.info(f"Generated information for {chemical}")
            return result.strip()
        
        except Exception as e:
            logger.error(f"Error getting chemical information: {str(e)}")
            return f"Unable to retrieve information for {chemical}: {str(e)}"
