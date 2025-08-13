# Chemical SMILES Chatbot

A Streamlit application that uses LangChain and a local LLaMA model to extract chemical names from natural language queries, fetch their SMILES representations from PubChem, and render their molecular structures.

## Features

- Natural language processing to extract chemical names
- SMILES (Simplified Molecular Input Line Entry System) representation lookup via PubChem
- Molecular structure rendering using RDKit
- Clean and responsive Streamlit UI

## Project Structure

```
chemical_smiles_chatbot/
├── app.py                  # Main Streamlit entry point
├── config/                 # Settings and constants
│   ├── __init__.py
│   ├── settings.py
│   └── prompts.py
├── src/                    # Source code
│   ├── __init__.py
│   ├── models/             # Code to load the local LLaMA model
│   │   ├── __init__.py
│   │   └── llama_model.py
│   ├── chains/             # Chains for chemical name extraction
│   │   ├── __init__.py
│   │   └── chemical_extractor.py
│   ├── services/           # SMILES lookup, molecule rendering
│   │   ├── __init__.py
│   │   ├── smiles_lookup.py
│   │   └── molecule_renderer.py
│   └── utils/              # Helpers like logging
│       ├── __init__.py
│       ├── logging_utils.py
│       └── smiles_utils.py
├── static/                 # CSS, JS, images
│   ├── css/
│   │   └── style.css
│   ├── js/
│   │   └── script.js
│   └── images/
├── templates/              # HTML templates
│   └── index.html
├── experiments/            # Jupyter Notebooks for testing
│   └── chemical_smiles_test.ipynb
├── data/                   # Local cache or sample data
│   └── smiles_cache.json
├── requirements.txt        # Dependencies
└── README.md               # This file
```

## Prerequisites

- Python 3.8+
- A local LLaMA model (e.g., llama-2-7b-chat.gguf)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/chemical_smiles_chatbot.git
cd chemical_smiles_chatbot
```

2. Install the required packages:
```bash
pip install -r requirements.txt
```

3. Download a LLaMA model and place it in the appropriate location (as specified in config/settings.py).

## Usage

1. Start the Streamlit application:
```bash
streamlit run app.py
```

2. Open your web browser and go to http://localhost:8501

3. Enter a chemical query in the text input field (e.g., "What is benzene?")

4. The application will:
   - Extract the chemical name from your query
   - Look up its SMILES representation
   - Render the molecular structure
   - Display the results

## Development

### Using a Mock Model for Development

If you don't have a local LLaMA model available, the application will use a mock model for development purposes. This allows you to test the application's functionality without the need for a large language model.

### Running the Experiments

You can run the Jupyter Notebook in the `experiments` folder to test the pipeline components individually:

```bash
jupyter notebook experiments/chemical_smiles_test.ipynb
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- PubChem for providing chemical data
- RDKit for molecular rendering capabilities
- LangChain for the LLM interaction framework
- Streamlit for the web application framework
