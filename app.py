import streamlit as st
import requests
import pandas as pd

st.set_page_config(page_title="Drug Solubility Predictor", layout="wide")

st.title("🧪 Pharma Drug Discovery Assistant")
st.subheader("Predicting Aqueous Solubility (logS)")

# API endpoint
API_URL = "http://127.0.0.1:8000/predict_solubility"

st.markdown("""
Enter a molecule's structure in **SMILES** format to predict its solubility.
Here are some examples:
- `CCO` (Ethanol)
- `c1ccccc1` (Benzene)
- `CC(=O)Oc1ccccc1C(=O)O` (Aspirin)
""")

smiles_input = st.text_input("Enter SMILES string:", "CC(=O)Oc1ccccc1C(=O)O")

if st.button("Predict Solubility"):
    if smiles_input:
        payload = {"smiles": smiles_input}
        try:
            response = requests.post(API_URL, json=payload)
            response.raise_for_status()  # Raise an exception for bad status codes
            result = response.json()
            st.success(f"**Predicted logS:** `{result['predicted_logS']:.4f}`")
        except requests.exceptions.RequestException as e:
            st.error(f"Error connecting to the API: {e}")
        except Exception as e:
            st.error(f"An error occurred: {e}")
    else:
        st.warning("Please enter a SMILES string.")