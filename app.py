import streamlit as st
import base64
import os
from io import BytesIO

# Import local modules
from src.models.llama_model import LlamaModel
from src.chains.chemical_extractor import ChemicalExtractor
from src.services.smiles_lookup import SMILESLookup
from src.services.molecule_renderer import MoleculeRenderer
from src.utils.logging_utils import setup_logger
from config.settings import APP_TITLE, APP_DESCRIPTION, MODEL_PATH

# Setup logging
logger = setup_logger(__name__)

def load_css():
    """Load custom CSS"""
    with open(os.path.join("static", "css", "style.css")) as f:
        st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)

def get_chemical_description(chemical_name):
    """
    Get a brief 1-2 sentence description of a chemical compound.
    
    Args:
        chemical_name (str): Name of the chemical
        
    Returns:
        str: Brief description of the chemical
    """
    chemical_descriptions = {
        "benzene": "Benzene is an organic compound with the formula C₆H₆, known for its distinctive ring structure. It's a colorless, highly flammable liquid with a sweet odor, widely used in industry as a precursor to other chemicals.",
        "toluene": "Toluene is an aromatic hydrocarbon with formula C₇H₈, consisting of a benzene ring with a methyl group. It's a colorless liquid commonly used as a solvent and in the production of other chemicals.",
        "acetaldehyde": "Acetaldehyde is an organic compound with formula C₂H₄O, produced naturally during metabolism. It's a colorless, volatile liquid with a pungent, fruity odor, used in various industrial processes.",
        "ethanol": "Ethanol is a simple alcohol with formula C₂H₆O, commonly known as drinking alcohol. It's a colorless, volatile liquid widely used as a solvent, fuel, and in beverages.",
        "methanol": "Methanol is the simplest alcohol with formula CH₄O, also known as wood alcohol. It's a colorless, toxic liquid used as a solvent, antifreeze, and fuel.",
        "methane": "Methane is the simplest hydrocarbon with formula CH₄, the main component of natural gas. It's a colorless, odorless gas that serves as a major energy source and greenhouse gas.",
        "water": "Water is a transparent, odorless, tasteless liquid with formula H₂O. It's essential for all known forms of life and covers about 71% of Earth's surface.",
        "carbon dioxide": "Carbon dioxide is a colorless gas with formula CO₂, produced by respiration and combustion. It's a crucial greenhouse gas and is used in photosynthesis by plants.",
        "acetone": "Acetone is an organic solvent with formula C₃H₆O, naturally produced in small amounts by the human body. It's a colorless, volatile liquid commonly used in nail polish remover and as an industrial solvent.",
        "caffeine": "Caffeine is a natural stimulant alkaloid found in coffee, tea, and many other plants. It's widely consumed to increase alertness and reduce fatigue.",
        "aspirin": "Aspirin is a medication used to reduce pain, fever, and inflammation. It's also used in low doses as a blood thinner to prevent heart attacks and strokes.",
        "glucose": "Glucose is a simple sugar with formula C₆H₁₂O₆, the primary source of energy for living organisms. It's found naturally in fruits and honey and is produced during photosynthesis.",
        "sodium chloride": "Sodium chloride, commonly known as table salt, has the formula NaCl. It's essential for life, used in food preparation, and plays important roles in biological processes.",
        "hydrogen peroxide": "Hydrogen peroxide is a chemical compound with formula H₂O₂, used as a disinfectant and bleaching agent. It naturally occurs in small amounts in the human body as a byproduct of metabolism.",
        "sulfuric acid": "Sulfuric acid is a highly corrosive strong mineral acid with formula H₂SO₄. It's widely used in industrial processes, including battery acid and chemical manufacturing.",
        "ammonia": "Ammonia is a compound with formula NH₃, consisting of nitrogen and hydrogen. It's a colorless gas with a distinctive pungent smell, commonly used in fertilizers and cleaning products."
    }
    
    # Return known description or generate a generic fallback
    return chemical_descriptions.get(chemical_name.lower(), 
                                   "This chemical is commonly used in research and industrial applications.")

def main():
    # Load the model and initialize components
    @st.cache_resource
    def load_model():
        logger.info("Loading LLaMA model...")
        model = LlamaModel()
        return model
    
    @st.cache_resource
    def init_components(_model):
        # Using leading underscore to tell Streamlit not to hash this argument
        extractor = ChemicalExtractor(_model)
        smiles_lookup = SMILESLookup()
        renderer = MoleculeRenderer()
        return extractor, smiles_lookup, renderer
    
    # Page configuration
    st.set_page_config(
        page_title=APP_TITLE,
        page_icon="🧪",
        layout="wide"
    )
    
    # Load custom CSS
    load_css()
    
    # App header
    st.title(APP_TITLE)
    st.markdown(APP_DESCRIPTION)
    
    # Initialize model and components
    with st.spinner("Loading model and components..."):
        model = load_model()
        extractor, smiles_lookup, renderer = init_components(model)
    
    # User input
    user_query = st.text_input("Enter your chemical query:", placeholder="What is benzene?")
    
    if user_query:
        with st.spinner("Processing your query..."):
            # Log the query for debugging
            logger.info(f"Processing query: {user_query}")
            
            try:
                # Extract chemical name from query
                chemical_name = extractor.extract_chemical_name(user_query)
                
                if chemical_name:
                    st.success(f"Detected chemical: {chemical_name}")
                    
                    # Show chemical description immediately
                    description = get_chemical_description(chemical_name)
                    st.info(f"**About {chemical_name}:** {description}")
                    
                    # Look up SMILES
                    try:
                        smiles = smiles_lookup.get_smiles(chemical_name)
                        st.info(f"SMILES: {smiles}")
                        
                        # Render molecule
                        if smiles:
                            img = renderer.render_molecule(smiles)
                            if img:
                                st.image(img, caption=f"Structure of {chemical_name}")
                    except Exception as e:
                        st.error(f"Error looking up chemical: {str(e)}")
                        logger.error(f"SMILES lookup error: {str(e)}")
                else:
                    st.warning("No specific chemical detected in your query. Please try again with a different query.")
                    logger.warning(f"No chemical detected in query: {user_query}")
            except Exception as e:
                st.error(f"Error processing query: {str(e)}")
                logger.error(f"Query processing error: {str(e)}", exc_info=True)
    
    # Additional information
    with st.expander("About this app"):
        st.markdown("""
        This app uses a local LLaMA model to extract chemical names from natural language queries,
        then fetches their SMILES representation from OPSIN, and renders the molecular structure.
        """)
        
    # Debug section
    if st.checkbox("Show debug information", False):
        st.subheader("Debug Information")
        st.write("This section is for troubleshooting the app.")
        
        with st.expander("Test chemical extraction"):
            test_query = st.text_input("Test query:", "What is benzene?")
            if st.button("Test extraction"):
                st.write(f"Testing extraction for: '{test_query}'")
                try:
                    # Use the LLM-based extraction
                    result = extractor.extract_chemical_name(test_query)
                    
                    if result:
                        st.success(f"Extracted chemical: {result}")
                    else:
                        st.warning("No chemical detected in the query")
                    
                except Exception as e:
                    st.error(f"Error during test: {str(e)}")

if __name__ == "__main__":
    main()
