# 🧪 Pharma Drug Discovery Assistant

A web application to predict the aqueous solubility (logS) of chemical compounds using a machine learning model. This tool is designed to assist in the early stages of drug discovery by providing quick solubility predictions from a molecule's SMILES string.

The project consists of two main components:
1.  A **FastAPI backend** that serves the pre-trained machine learning model.
2.  A **Streamlit frontend** that provides a user-friendly web interface to interact with the model.

 <!-- It's a good idea to add a screenshot of your app here -->

---

## 🚀 Features

-   **Solubility Prediction**: Predicts the logarithm of the aqueous solubility (logS) for a given molecule.
-   **SMILES Input**: Accepts molecule structures in the widely-used SMILES format.
-   **Interactive UI**: A simple and clean web interface built with Streamlit.
-   **REST API**: A robust API endpoint for programmatic access to the prediction model.
-   **Fast & Efficient**: Built with FastAPI for high performance.

---

## 🏛️ Project Architecture

The application is decoupled into a backend service and a frontend client:

-   **Backend (`main.py`)**:
    -   Built with **FastAPI**.
    -   Loads a pre-trained `scikit-learn` model (`solubility_model.joblib`).
    -   Uses **RDKit** to convert SMILES strings into Morgan fingerprints (molecular features).
    -   Exposes a `/predict_solubility` endpoint to make predictions.

-   **Frontend (`app.py`)**:
    -   Built with **Streamlit**.
    -   Provides a text input for users to enter a SMILES string.
    -   Communicates with the FastAPI backend via HTTP requests to get predictions.
    -   Displays the predicted solubility to the user.

---

## 🛠️ Tech Stack

-   **Backend**: Python, FastAPI, Uvicorn
-   **Frontend**: Streamlit
-   **ML & Data**: Scikit-learn, RDKit, Pandas, NumPy
-   **Model Training Data**: Delaney (ESOL) Dataset

---

## ⚙️ Installation & Setup

Follow these steps to set up and run the project locally.

### 1. Clone the Repository

```bash
git clone https://github.com/sudheersivakumar/Pharma-Drug-Discovery-Assistant
cd pharma-drug-discovery-assist
```

### 2. Create a Virtual Environment (Recommended)

```bash
# For Windows
python -m venv venv
venv\Scripts\activate

# For macOS/Linux
python3 -m venv venv
source venv/bin/activate
```

### 3. Install Dependencies

All required packages are listed in `requirements.txt`.

```bash
pip install -r requirements.txt
```

---

## ▶️ Running the Application

You need to run the backend and frontend servers in two separate terminals.

### 1. Start the FastAPI Backend

In your first terminal, run the following command from the project root directory:

```bash
uvicorn main:app --reload
```

The API will be available at `http://127.0.0.1:8000`. You can access the interactive API documentation at `http://127.0.0.1:8000/docs`.

### 2. Start the Streamlit Frontend

In your second terminal, run this command:

```bash
streamlit run app.py
```

The web application will open in your browser, and you can access it at `http://localhost:8501`.

---

## 📂 Project Structure

```
.
├── app.py                  # Streamlit frontend application
├── main.py                 # FastAPI backend API
├── requirements.txt        # Python dependencies
├── solubility_model.joblib # Pre-trained solubility prediction model
├── README.md               # This file
└── archive/
    └── delaney.csv         # Dataset used for training the model
```

---

## 🔌 API Usage

You can also interact with the API directly.

### Endpoint: `/predict_solubility`

-   **Method**: `POST`
-   **Request Body**: A JSON object with the SMILES string.
    ```json
    {
      "smiles": "CC(=O)Oc1ccccc1C(=O)O"
    }
    ```
-   **Example with `curl`**:
    ```bash
    curl -X POST "http://127.0.0.1:8000/predict_solubility" \
    -H "Content-Type: application/json" \
    -d '{"smiles": "CC(=O)Oc1ccccc1C(=O)O"}'
    ```
-   **Success Response (200 OK)**:
    ```json
    {
      "smiles": "CC(=O)Oc1ccccc1C(=O)O",
      "predicted_logS": -3.66  // Example value
    }
    ```

---

## 📄 License

This project is licensed under the MIT License. See the `LICENSE` file for details.

