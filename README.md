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
- 
