from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
import joblib
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# Initialize FastAPI app
app = FastAPI(
    title="Pharma Drug Discovery Assistant API",
    description="An API to predict molecular properties.",
    version="0.1.0",
)

# Load the trained model
model = joblib.load('C:/Users/sudhe/OneDrive/Documents/Local RAG/Pharma drug discovery Assist/solubility_model.joblib')

class Molecule(BaseModel):
    smiles: str

def smiles_to_fingerprint(smiles: str):
    """Converts a SMILES string to a Morgan fingerprint."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    return np.array(fp)

@app.post("/predict_solubility")
def predict_solubility(molecule: Molecule):
    """
    Predicts the water solubility (logS) of a molecule from its SMILES string.
    """
    fingerprint = smiles_to_fingerprint(molecule.smiles)
    if fingerprint is None:
        raise HTTPException(status_code=400, detail="Invalid SMILES string provided.")

    # Model expects a 2D array
    prediction = model.predict(fingerprint.reshape(1, -1))

    return {"smiles": molecule.smiles, "predicted_logS": prediction[0]}

@app.get("/")
def read_root():
    return {"message": "Welcome to the Drug Discovery Assistant API"}