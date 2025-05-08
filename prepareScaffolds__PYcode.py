import argparse
import base64
import csv
import io
import json
import re
import signal
import sys
import time
import zlib
from collections import Counter
from func_timeout import func_timeout, FunctionTimedOut
from operator import itemgetter
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import Geometry

## Calculate scaffolds
	#	1.		Get the structures in SMILES from TAB
	#	1.1.		Get mol
	#	1.2.		Get Canonical Smiles
	#	1.3			Make unique
	#######
	#	2.		Get mol
	#	3.		Get number of components
	#	4.		Get the largest component
	#	5.		Get the MW
	#	6.		Check MW (1800 for now)
	#	7.		Get number of CC bonds
	#	8.		Check number of CC bonds (1 for now)
	#	9.		Remove stereo, Hs
	#	10.		Get scaff
	#	11.		Repeating SMILES: Count the number of SMILES repeating other SMILES: if there are 2 exact matches -> 1 repeat; 3 exact matches -> 2 repeats, 2 & 3 distinct exact matches -> 3 repeats
	##Some vars to store the intermediate results
initial_list = []
result_list  = []
with open(".../input/chembl_pubchem_25-04-25.tsv", newline='') as tab_file:
	smiles_reader = csv.reader(tab_file, delimiter='	')
	for row in smiles_reader:
		current_mol = 'NONE'
		original_smiles = row[1]
		try:
			current_mol = Chem.MolFromSmiles(original_smiles)
		except Exception:
			continue
		if current_mol == 'NONE':
			continue
		try:
			smiles = Chem.MolToSmiles(current_mol)
		except Exception:
			continue
		initial_list.append(smiles)

initial_set = set(initial_list)
smiles_list = list(initial_set)
for record in smiles_list:
	current_mol = 'NONE'
	current_mol = Chem.MolFromSmiles(record)
	if current_mol == 'NONE':
		continue
	n_components = len(Chem.rdmolops.GetMolFrags(current_mol))
	if n_components > 1:
		current_mol = rdMolStandardize.LargestFragmentChooser().choose(current_mol)
	if Descriptors.MolWt(current_mol) > 1800:
		continue
	n_cc_bonds = 0
	n_bonds = 0
	for bond in current_mol.GetBonds():
			n_bonds += 1
			anum1 = bond.GetBeginAtom().GetAtomicNum()
			anum2 = bond.GetEndAtom().GetAtomicNum()
			if anum1 == 6 and anum2 == 6:
				n_cc_bonds += 1
	if n_cc_bonds < 1:
		continue
	scaff = ''
	scaff_gen = ''
	scaff_smiles = 'NONE'
	scaff_gen_smiles = 'NONE'
	dummy_plain = Chem.RemoveStereochemistry(current_mol)
	current_mol = Chem.RemoveHs(current_mol)
	smiles = Chem.MolToSmiles(current_mol)
	try:
		scaff = MurckoScaffold.GetScaffoldForMol(current_mol)
		scaff_gen = MurckoScaffold.GetScaffoldForMol(MurckoScaffold.MakeScaffoldGeneric(scaff))
		scaff_smiles = Chem.MolToSmiles(scaff)
		scaff_gen_smiles = Chem.MolToSmiles(scaff_gen)
	except Exception:
		continue
	if scaff_smiles == '': scaff_smiles = '*'
	if scaff_smiles == 'NONE': scaff_smiles = '[*-]'
	if scaff_gen_smiles == '': scaff_gen_smiles = '*'
	if scaff_gen_smiles == 'NONE': scaff_gen_smiles = '[*-]'
	result = [record, scaff_smiles, scaff_gen_smiles]
	result_list.append(result)



with open('.../output/chembl_pubchem_25-04-25_scaff_double.tsv', 'w', newline='') as file:
	writer = csv.writer(file, delimiter='\t')
	writer.writerows(result_list)



