from rdkit import Chem
import torch
from torch_geometric.data import Data
import pandas as pd

def convert_smiles_to_graph(string):
    
    molecule = Chem.MolFromSmiles(string)
    if molecule is None:
        print("Molecule not found")
        return None

    #Nodes = atoms
    nodes = []

    #Gets the atomic number for each element in the molecule 
    #not just the letter (since it can't be interpreted well by the
    #computer) 
    for atom in molecule.GetAtoms():
        nodes.append([atom.GetAtomicNum()])

    #Converts the node array into a pytorch tensor (which is pretty much an array, but
    #multi-dimensional and for torch specifically) 
    x=torch.tensor(nodes, dtype = torch.float) 

    #Edges = bonds 

    #Contains all indeces of atoms, from which the bond starts 
    edge_start = []
    #Contains all indeces of atoms, from which the bond starts 
    edge_end = []
    #Contains the bond types for each bond
    edge_types = []

    for bond in molecule.GetBonds():
        first_atom = bond.GetBeginAtomIdx()
        second_atom = bond.GetEndAtomIdx()

        #1.0 = single, 2.0 = double, 3.0 = triple, 1.5 = aromatic
        b_type = bond.GetBondTypeAsDouble()

        #Appends twice for the 2 variations of the bond

        edge_types.append([b_type])
        edge_types.append([b_type])
        
        #Adds both combinations, since bonds don't have directions 
        edge_start.extend([first_atom, second_atom])
        edge_end.extend([second_atom, first_atom])

    #Matrix connecting edge starts and ends
    edge_index = torch.tensor([edge_start, edge_end], dtype = torch.long)
    edge_attr = torch.tensor(edge_types, dtype = torch.float)

    graph_data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr)
    return graph_data

df = pd.read_csv("antidepressant_smiles.csv")
drug_graphs = {}

#Traverse through each row
for i, row in df.iterrows():
    drug_name = row["Drug:"]
    smiles = row["SMILES"]

    graph = convert_smiles_to_graph(smiles)

    if graph is not None:
        drug_graphs[drug_name] = graph

    else:
        print(f"Graph not done for {drug_name}")

print("Complete")
