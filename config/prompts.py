"""
Define prompts used by the LangChain components.
"""

CHEMICAL_EXTRACTION_TEMPLATE = """
You are a professional chemist specializing in chemical identification. 

User query: {query}

Your task: Extract ONLY the chemical name from the query above.

Examples:
"What is benzene?" → benzene
"Tell me about aspirin" → aspirin
"What are the properties of H2O?" → water
"Can you explain how caffeine works?" → caffeine
"What's the chemical formula for table salt?" → sodium chloride
"What is carbon di oxide?" → carbon dioxide
"Tell me about Hydrogen-Peroxide" → hydrogen peroxide
"What is CH3COOH?" → acetic acid
"What is benzaldehyde" → benzaldehyde

Rules:
1. Return ONLY the chemical name, nothing else
2. Standardize the name (e.g., "carbon di oxide" → "carbon dioxide")
3. Convert formulas to names when common (H2O → water)
4. If no chemical is mentioned, return "None"

Chemical name:
"""

CHEMICAL_INFORMATION_TEMPLATE = """
You are a chemistry expert assistant. Provide information about the following chemical:

Chemical: {chemical}
SMILES: {smiles}

Provide a brief description of this chemical, including:
1. Basic properties (molecular formula, molecular weight, state at room temperature)
2. Common uses and applications
3. Important characteristics and reactivity
4. Safety information and common hazards (if applicable)
5. Interesting facts or historical significance

Format your answer as a concise but informative paragraph that would be helpful to someone learning about this chemical.

Information:
"""
