import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
import joblib
import numpy as np
from joblib import Parallel, delayed


def _get_fingerprint(mol):
    """Helper function to calculate a single fingerprint. This is needed to make the call picklable for joblib."""
    # Use the modern MorganGenerator API to avoid deprecation warnings
    fp_gen = AllChem.GetMorganGenerator(radius=2, fpSize=2048)
    return fp_gen.GetFingerprintAsNumPy(mol)


def generate_fingerprints(mols, n_jobs=-1):
    """Generate Morgan fingerprints for a list of mol objects in parallel."""
    tasks = (delayed(_get_fingerprint)(mol) for mol in mols)
    return Parallel(n_jobs=n_jobs, verbose=1)(tasks)


def train_model():
    """Loads data, trains a model, and saves it."""
    print("Loading data...")
    df = pd.read_csv('C:/Users/sudhe/OneDrive/Documents/Local RAG/Pharma drug discovery Assist/archive/delaney.csv')

    # The standard Delaney dataset uses 'measured log solubility in mols per litre'
    # for the target variable. We'll rename it to 'logS' for convenience.
    # If your column has a different name, please update it here.
    solubility_col_original = 'measured log(solubility:mol/L)'
    if solubility_col_original in df.columns:
        df.rename(columns={solubility_col_original: 'logS'}, inplace=True)

    print(f"Available columns: {df.columns.tolist()}")

    # Create mol objects from SMILES and filter out invalid ones
    df['mol'] = df['SMILES'].apply(Chem.MolFromSmiles)
    original_count = len(df)
    df.dropna(subset=['mol'], inplace=True)
    if len(df) < original_count:
        print(f"Warning: Removed {original_count - len(df)} invalid SMILES strings.")

    print("Generating fingerprints...")
    # Generate fingerprints in parallel
    X = np.array(generate_fingerprints(df['mol'].tolist()))
    y = df['logS'].values

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    print("Training model...")
    # Use n_jobs=-1 to use all available CPU cores for training
    model = RandomForestRegressor(n_estimators=100, random_state=42, n_jobs=-1)
    model.fit(X_train, y_train)

    print("Evaluating model...")
    y_pred = model.predict(X_test)
    mse = mean_squared_error(y_test, y_pred)
    print(f"Model MSE: {mse:.4f}")

    print("Saving model...")
    joblib.dump(model, 'solubility_model.joblib')
    print("Model saved as solubility_model.joblib")

if __name__ == "__main__":
    train_model()
