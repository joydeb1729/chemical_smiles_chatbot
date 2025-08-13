"""
Module for rendering molecule structures from SMILES.
"""

import io
import base64
try:
    from rdkit import Chem
    from rdkit.Chem import Draw, AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    from PIL import Image, ImageDraw, ImageFont

from src.utils.logging_utils import setup_logger

logger = setup_logger(__name__)

class MoleculeRenderer:
    """
    Class to render molecular structures from SMILES representations using RDKit.
    """
    
    def __init__(self, image_size=(400, 300), use_svg=False):
        """
        Initialize the molecule renderer.
        
        Args:
            image_size (tuple, optional): Size of the rendered image. Defaults to (400, 300).
            use_svg (bool, optional): Whether to use SVG format. Defaults to False.
        """
        self.image_size = image_size
        self.use_svg = use_svg
        self.rdkit_available = RDKIT_AVAILABLE
        if not self.rdkit_available:
            logger.warning("RDKit is not available. Will use fallback rendering.")
    
    def render_molecule(self, smiles, highlight_atoms=None):
        """
        Render a molecule from its SMILES representation.
        
        Args:
            smiles (str): The SMILES representation of the molecule.
            highlight_atoms (list, optional): List of atom indices to highlight. Defaults to None.
            
        Returns:
            bytes: The rendered image as bytes, or None if rendering failed.
        """
        if not smiles:
            logger.error("Cannot render molecule: No SMILES provided")
            return None
            
        try:
            if self.rdkit_available:
                return self._render_with_rdkit(smiles, highlight_atoms)
            else:
                return self._render_fallback(smiles)
        except Exception as e:
            logger.error(f"Error rendering molecule: {str(e)}")
            return self._render_fallback(smiles, error=str(e))
    
    def _render_with_rdkit(self, smiles, highlight_atoms=None):
        """
        Render a molecule using RDKit.
        
        Args:
            smiles (str): The SMILES representation of the molecule.
            highlight_atoms (list, optional): List of atom indices to highlight.
            
        Returns:
            Image: The rendered image.
        """
        logger.info(f"Rendering molecule for SMILES: {smiles}")
        
        # Parse the SMILES string to create a molecule object
        mol = Chem.MolFromSmiles(smiles)
        
        if not mol:
            logger.warning(f"Failed to parse SMILES: {smiles}")
            return self._render_fallback(smiles, error="Invalid SMILES structure")
        
        # Add hydrogens to the molecule for a more complete representation
        mol = Chem.AddHs(mol)
        
        # Generate 2D coordinates for the molecule
        AllChem.Compute2DCoords(mol)
        
        # Render the molecule
        if self.use_svg:
            # SVG rendering
            svg = Draw.MolToImage(mol, size=self.image_size, highlightAtoms=highlight_atoms)
            img_data = io.BytesIO()
            svg.save(img_data, format="SVG")
            return img_data.getvalue()
        else:
            # PNG rendering
            img = Draw.MolToImage(mol, size=self.image_size, highlightAtoms=highlight_atoms)
            return img
    
    def _render_fallback(self, smiles, error=None):
        """
        Fallback rendering when RDKit is not available or fails.
        
        Args:
            smiles (str): The SMILES representation of the molecule.
            error (str, optional): Error message to display. Defaults to None.
            
        Returns:
            Image: A simple image with the SMILES text.
        """
        # Create a placeholder image with the SMILES string
        img = Image.new('RGB', self.image_size, color='white')
        d = ImageDraw.Draw(img)
        
        # Add the SMILES string
        d.text((10, 10), f"SMILES: {smiles}", fill='black')
        
        # Add error message if provided
        if error:
            d.text((10, 50), f"Error: {error}", fill='red')
            d.text((10, 90), "Note: RDKit is required for proper molecular rendering.", fill='blue')
        else:
            d.text((10, 50), "RDKit is required for proper molecular rendering.", fill='blue')
        
        return img
    
    def get_molecule_as_base64(self, smiles):
        """
        Get a molecule image as a base64-encoded string for embedding in HTML.
        
        Args:
            smiles (str): The SMILES representation of the molecule.
            
        Returns:
            str: Base64-encoded image string, or None if rendering failed.
        """
        try:
            img = self.render_molecule(smiles)
            
            if not img:
                return None
            
            # Convert the image to base64
            if isinstance(img, Image.Image):
                buffered = io.BytesIO()
                img.save(buffered, format="PNG")
                img_str = base64.b64encode(buffered.getvalue()).decode()
                return f"data:image/png;base64,{img_str}"
            else:
                # Already bytes (SVG)
                img_str = base64.b64encode(img).decode()
                return f"data:image/svg+xml;base64,{img_str}"
        except Exception as e:
            logger.error(f"Error converting molecule to base64: {str(e)}")
            return None
            
        except Exception as e:
            logger.error(f"Error converting molecule to base64: {str(e)}")
            return None
