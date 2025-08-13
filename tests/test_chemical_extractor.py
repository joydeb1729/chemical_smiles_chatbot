import pytest
from src.chains.chemical_extractor import ChemicalExtractor
from unittest.mock import MagicMock

class TestChemicalExtractor:
    def setup_method(self):
        # Create a mock LLM
        self.mock_llm = MagicMock()
        self.mock_llm.model = MagicMock()
        self.mock_llm.model.return_value = "test_chemical"
        
        # Create the extractor with the mock LLM
        self.extractor = ChemicalExtractor(self.mock_llm)
    
    def test_extract_benzene(self):
        # Test direct extraction without LLM for "benzene"
        result = self.extractor.extract_chemical_name("What is benzene?")
        assert result.lower() == "benzene"
        
        result = self.extractor.extract_chemical_name("Tell me about benzene")
        assert result.lower() == "benzene"
    
    def test_extract_water(self):
        # Test extraction for water with different forms
        result = self.extractor.extract_chemical_name("What is H2O?")
        assert result.lower() == "water"
        
        result = self.extractor.extract_chemical_name("Properties of water")
        assert result.lower() == "water"
    
    def test_extract_with_llm(self):
        # Test extraction that requires the LLM
        result = self.extractor.extract_chemical_name("What is the structure of dimethyl sulfoxide?")
        # This should use the mock LLM which returns "test_chemical"
        assert result == "test_chemical"
    
    def test_no_chemical(self):
        # Mock the LLM to return "None" for this test
        self.mock_llm.model.return_value = "None"
        
        # Test query with no chemical
        result = self.extractor.extract_chemical_name("What is the weather like today?")
        assert result is None
