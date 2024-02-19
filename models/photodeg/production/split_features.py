#!/usr/bin/python3

import pandas as pd
import sys

input_pickle = sys.argv[1]
atom_pickle = sys.argv[2]
bond_pickle = sys.argv[3]

all_features = pd.read_pickle(input_pickle)
smiles_indexed = all_features.set_index('smiles')
atom_features = smiles_indexed[['partial_charge', 'partial_neu', 'partial_elec', 'NMR']]
bond_features = smiles_indexed[['bond_order', 'bond_distance']]

atom_features.to_pickle(atom_pickle)
bond_features.to_pickle(bond_pickle)