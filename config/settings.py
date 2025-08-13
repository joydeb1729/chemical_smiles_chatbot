"""
Configuration settings for the Chemical SMILES Chatbot application.
"""

# Application settings
APP_TITLE = "Chemical SMILES Chatbot"
APP_DESCRIPTION = "Ask questions about chemicals and get their SMILES representation and structure."

# Model settings
MODEL_PATH = "models/llama-2-7b-chat.Q4_K_M.gguf"  # Path to your local LLaMA model
MODEL_TEMP = 0.1
MODEL_MAX_TOKENS = 256
MODEL_N_CTX = 2048  # Context window size
MODEL_N_GPU_LAYERS = 0  # Number of layers to offload to GPU (0 for CPU only)

# Service settings
PUBCHEM_TIMEOUT = 10  # Timeout for PubChem API calls in seconds
CACHE_EXPIRY = 86400  # Cache expiry time in seconds (24 hours)

# Logging settings
LOG_LEVEL = "INFO"
LOG_FILE = "app.log"

# Path settings
DATA_DIR = "data"
CACHE_DIR = f"{DATA_DIR}/cache"
