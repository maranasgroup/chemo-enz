### Requirements: 

1. Rdkit
2. Streamlit
3. Pandas
4. Numpy
5. Keras
6. matplotlib
7. Pulp
8. CPLEX solver
    

### The files associated with the MinChemBio are seperately attached here, since the github limits uploading files larger than 100 MB without LFS.
[Dataset and code on Scholarsphere](https://scholarsphere.psu.edu/resources/28575451-cb5a-47db-a9e9-d42c984b9ddc)


# How to use minChemBio tool


## Getting Started

1. **Create the Conda environment**  
   Set up the environment with all required dependencies by running:
   ```bash
   conda env create -f minchembio_yml.yml
   ```

2. **Prepare the dataset**  
   Ensure you have the dataset files in the working directory

3. **Extract reactant and product IDs**  
   Open the provided Jupyter notebook and run it to extract the IDs of reactants and products from `S1.txt`.

4. **Run the Streamlit web app**  
   Launch the web app interface by executing the following command:
   ```bash
      streamlit run minchembio_streamlit.py
   ```
   OR

    ```bash
      python run minchembio.py
   ```
   Enter the reactant and product IDs in the app to generate solutions using **minChemBio**.

6. **Visualize the solutions**  
   Use the visualization Jupyter notebook to explore the generated pathways. Provide the folder containing the pathway outputs as input to the notebook.

7. **Reference data (optional)**  
   The file `S2.csv` contains a list of USPTO patent IDs and their corresponding reaction SMILES strings. This can be used for further analysis or comparison.



## File List

This repository contains the code and datasets required to run **minChemBio**, a tool for exploring pathways using mixed-integer linear programming (MILP).

---

### Data Files

1. **`all_rij_with_miss_cat.json`**  
   A dictionary where **molecule IDs** are the keys. Each value is a dictionary listing all reactions involving that molecule, either as a **reactant (-1)** or a **product (1)**.

2. **`all_sij_with_miss_cat.json`**  
   A dictionary where **reaction IDs** are the keys. Each value is a dictionary listing all molecules involved in the reaction as **reactants (-1)** or **products (1)**.

3. **`bio_chem_smiles_ids_dict_NEW.json`**  
   A dictionary with **canonical SMILES strings** as keys and **molecule IDs** as values.

4. **`bio_chem_smiles_ids_dict_updated.json`**  
   A dictionary with **canonical SMILES strings** as keys and **molecule IDs** as values. Similar to earlier SMILES-ID mappings but includes fewer molecules.

5. **`bio_chem_ids_dict_NEW.json`**  
   A dictionary with **molecule IDs** as keys and **canonical SMILES strings** as values. Similar to earlier SMILES-ID mappings but includes fewer molecules.

6. **`metanetx_metab_db.json`** 
   A dictionary where **MetaNetx molecule  IDs** are the keys. Each value is a dictionary listing _Name_, _Formula_, _Charge_, _Mass_, _InChI_ , _InchIKey_ , _SMILES_ , _Reference_ for each molecule.
   This dataset is extracted from the [MetaNetX web platform](https://www.metanetx.org/mnxdoc/mnxref.html).

7. **`rev_pair_90_nondup.json`**  
   A dictionary where **reaction IDs** are the keys. Each value is a list containing reverse mappings extracted from the same reaction. Used to ensure that **forward and reverse reactions do not co-occur** in the same pathway.

8. **`rxn_classify_with_miss_cat.json`**  
   A dictionary with **reaction IDs** as keys. The value for each key is:
   - `1` for **chemical reactions**
   - `0` for **biological reactions**

9. **`S1.txt`**  
   A text file containing **molecule IDs** and their corresponding **canonical SMILES strings**.

10. **`S2.csv`**  
   A CSV file containing **reaction IDs** from the USPTO dataset, along with their associated **patent number, year**, and **reaction SMILES string**.

---

### Code Files

11. **`minchembio_streamlit.py`**  
   A Streamlit web app interface for **minChemBio**.  
   - **Inputs**: Product and reactant molecule IDs  
   - **Required files**:  
     - `all_rij_with_miss_cat.json`  
     - `all_sij_with_miss_cat.json`  
     - `bio_chem_smiles_ids_dict_updated.json`  
     - `rev_pair_90_nondup.json`  
     - `rxn_classify_with_miss_cat.json`  
   - **Output**: A text file named in the format  
     `productID_from_reactionID-timestamp_.txt`  
     This file contains all possible solutions (pathways), each being a list of **reaction IDs** derived from solving the MILP problem.

12. **`minchembio.py`**  
   A Python script version of the Streamlit app.  
   - Same functionality as `minchembio_streamlit.py`  
   - Users need to **edit the `main()` function** to input the desired molecule IDs.

13. **`visualize.ipynb`**  
   A Jupyter notebook for visualizing the output pathways.  
   - **Inputs**: Same as `minchembio_streamlit.py` + the results file generated from MILP  
   - **Output**: Visual representations of all identified pathways, saved as **.png** files.



