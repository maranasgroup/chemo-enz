### Requirements: 

1. Rdkit
2. Streamlit
3. Pandas
4. Numpy
5. Keras
6. matplotlib
7. Pulp
8. CPLEX solver
    

The files associated with the MinChemBio are seperately attached here, since the github limits uploading files larger than 100 MB without LFS.
[Dataset and code on Scholarsphere](https://doi.org/10.26207/tbg0-gr88)


# How to use minChemBio tool

```python
conda env create -f environment.yaml
```
- Get the dataset ready (S1.txt)
- Now get the ids of reactant and product using the jupyter code
- Use the ids to put into the streamlit web app by executing the file minchembio_streamlit.py
- Now, you have the solutions given by minChemBio
- To visualize these solutions, you have to use the another juptyer code, by inputing the folder of the pathways
- S2.csv is the list of patent id of USPTO reactions with their reaction SMILES strings

Sure! Here's the updated section in Markdown syntax, ready to paste directly into your `README.md` file:

```markdown
## Getting Started

1. **Create the Conda environment**  
   Set up the environment with all required dependencies by running:
   ```bash
   conda env create -f environment.yaml
   ```

2. **Prepare the dataset**  
   Ensure you have the dataset file named `S1.txt` in the working directory.

3. **Extract reactant and product IDs**  
   Open the provided Jupyter notebook and run it to extract the IDs of reactants and products from `S1.txt`.

4. **Run the Streamlit web app**  
   Launch the web app interface by executing the following command:
   ```bash
   python minchembio_streamlit.py
   ```
   Enter the reactant and product IDs in the app to generate solutions using **minChemBio**.

5. **Visualize the solutions**  
   Use the visualization Jupyter notebook to explore the generated pathways. Provide the folder containing the pathway outputs as input to the notebook.

6. **Reference data (optional)**  
   The file `S2.csv` contains a list of USPTO patent IDs and their corresponding reaction SMILES strings. This can be used for further analysis or comparison.
```

Let me know if you want me to add a **Project Overview**, **Installation**, or **License** section too!
