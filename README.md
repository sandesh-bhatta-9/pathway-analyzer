# Interactive Pathway Overlap Analyzer

This tool is a web application for analyzing genetic overlaps between biological pathways from the KEGG database. It is designed to help researchers identify shared molecular mechanisms between various cancers, metabolic diseases like diabetes, and drug metabolism pathways, making it suitable for chemical and cancer informatics research.

**Live App Demo:** [https://pathway-analyzer-y3oprqvkoezrqdydnjlps5.streamlit.app](https://pathway-analyzer-y3oprqvkoezrqdydnjlps5.streamlit.app)

## Features

-   **Automated Pathway Fetching**: Connects directly to the KEGG API to pull the most up-to-date list of all human pathways.
-   **Interactive Selection**: Users can select multiple pathways of interest from a filterable list in the sidebar.
-   **Overlap Analysis**: Generates a table of shared genes, ranked by how frequently they appear across the selected pathways.
-   **Network Visualization**: Creates a dynamic, color-coded network graph to illustrate the connections between different pathway types and their shared genes.



## How to Use the Live App

1.  Navigate to the [live app link](https://pathway-analyzer-y3oprqvkoezrqdydnjlps5.streamlit.app).
2.  Use the multiselect box in the sidebar to choose the pathways you wish to analyze.
3.  The results, including the shared gene table and network graph, will update automatically on the main screen.

---

## Running the Code Locally

To run this application on your own machine, follow these steps:

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/sandesh-bhatta-9/pathway-analyzer.git](https://github.com/sandesh-bhatta-9/pathway-analyzer.git)
    cd pathway-analyzer
    ```

2.  **Install the required libraries:**
    ```bash
    pip install -r requirements.txt
    ```

3.  **Run the Streamlit app:**
    ```bash
    streamlit run app.py
    ```
