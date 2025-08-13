"""
Define prompts used by the LangChain components.
"""

CHEMICAL_EXTRACTION_TEMPLATE = """
You are a helpful assistant that extracts chemical names from text.

User query: {query}

Your task is to identify and extract the name of a chemical compound from the query.
Common query patterns include "What is [chemical]?", "Tell me about [chemical]", etc.

Examples:
- Query: "What is benzene?"
  Output: benzene
- Query: "Tell me about aspirin"
  Output: aspirin
- Query: "What are the properties of H2O?"
  Output: water
- Query: "Can you explain how caffeine works?"
  Output: caffeine
- Query: "What's the chemical formula for table salt?"
  Output: sodium chloride
- Query: "What is carbon di oxide?"
  Output: carbon dioxide
- Query: "Tell me about Hydrogen-Peroxide"
  Output: hydrogen peroxide

If there are multiple chemicals, extract the main one that is the focus of the query.
If there is no specific chemical mentioned, return None.

Important:
1. Standardize chemical names (e.g., "carbon di oxide" → "carbon dioxide", "table salt" → "sodium chloride")
2. Return ONLY the standardized chemical name without any explanations
3. Convert chemical formulas to names when possible (e.g., H2O → water, CO2 → carbon dioxide)
4. Correct common misspellings and spacing variations

Chemical name:
"""

CHEMICAL_INFORMATION_TEMPLATE = """
You are a chemistry expert assistant. Provide information about the following chemical:

Chemical: {chemical}
SMILES: {smiles}

Provide a brief description of this chemical, including:
1. Basic properties
2. Common uses
3. Important characteristics

Information:
"""
