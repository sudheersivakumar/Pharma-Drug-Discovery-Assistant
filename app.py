import streamlit as st
import requests
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
from PIL import Image
import io

def smiles_to_image(smiles: str):
    """Generates a 2D image of a molecule from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Draw.MolToImage(mol, size=(400, 400))

def calculate_properties(smiles: str):
    """Calculates key physicochemical properties from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    return mw, logp

def interpret_logS(logS: float):
    """Provides a qualitative interpretation of the logS value."""
    if logS > -2:
        return "High", "🟢"
    elif -4 <= logS <= -2:
        return "Moderate", "🟡"
    else:
        return "Low", "🔴"

st.set_page_config(page_title="Drug Solubility Predictor", layout="wide")

# Custom CSS for a more modern and clean look
st.markdown("""
<style>
    .main {
        background-color: #F8F9FA;
    }
    .st-emotion-cache-16txtl3 {
        padding: 2rem 2rem 10rem;
    }
    .stButton>button {
        border-radius: 8px;
        padding: 10px 24px;
        font-weight: bold;
    }
    .stTextInput>div>div>input {
        border-radius: 8px;
    }
    div[data-testid="stMetric"] {
        background-color: #000000;
        border: 1px solid #E0F7FA;
        border-radius: 8px;
        padding: 15px;
    }
</style>
""", unsafe_allow_html=True)

st.title("🧪 Pharma Drug Discovery Assistant")
st.markdown("### A Modern Tool for Predicting Molecular Properties")

# API endpoint
API_URL = "http://127.0.0.1:8000/predict_solubility"

col1, col2 = st.columns([1, 1.5])

with col1:
    st.header("Molecule Input")
    st.markdown("Enter a molecule's structure in **SMILES** format.")

    # Example molecules
    example_smiles = {
        "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
        "Caffeine": "Cn1cnc2n(C)c(=O)n(C)c(=O)c12",
        "Ibuprofen": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
        "Paracetamol": "CC(=O)Nc1ccc(O)cc1",
    }
    selected_example = st.selectbox("Or choose an example:", list(example_smiles.keys()))
    smiles_input = st.text_input("Enter SMILES string:", example_smiles[selected_example])

    predict_button = st.button("✨ Predict Properties")

with col2:
    st.header("Prediction Results")
    if predict_button and smiles_input:
        with st.spinner('Analyzing molecule and predicting properties...'):
            payload = {"smiles": smiles_input}
            try:
                # API call for solubility
                response = requests.post(API_URL, json=payload)
                response.raise_for_status()
                result = response.json()
                predicted_logS = result['predicted_logS']

                # Calculate other properties and get image
                mw, logp = calculate_properties(smiles_input)
                mol_image = smiles_to_image(smiles_input)

                if mol_image and mw is not None:
                    st.subheader("Predicted Aqueous Solubility (logS)")
                    interpretation, emoji = interpret_logS(predicted_logS)
                    st.metric(label=f"Solubility Level: {interpretation} {emoji}", value=f"{predicted_logS:.4f}")

                    st.subheader("Molecule Structure")
                    st.image(mol_image, caption=f"2D Structure for: {smiles_input}", use_column_width=True)

                    st.subheader("Other Physicochemical Properties")
                    prop_data = {"Property": ["Molecular Weight (g/mol)", "LogP (Lipophilicity)"], "Value": [f"{mw:.2f}", f"{logp:.2f}"]}
                    st.table(pd.DataFrame(prop_data))
                else:
                    st.error("Invalid SMILES string. Could not generate molecule structure or properties.")
            except requests.exceptions.RequestException as e:
                st.error(f"API Error: Could not connect to the prediction service. Please ensure the backend is running.")
            except Exception as e:
                st.error(f"An unexpected error occurred: {e}")
    else:
        st.info("Enter a SMILES string and click 'Predict Properties' to see the results here.")
