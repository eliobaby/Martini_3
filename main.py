from rdkit import Chem
from setup import get_atom_properties, connectivity_matrix, duplicate_ring_number_handler, remove_stars
from martini_3_dictionary import get_m3_dict
from mapping_scheme import map_molecule, parse_smiles
from algorithm_c1 import map_martini_beads
from text_file import gro_maker, itp_maker, export_bead_mapping
import copy

# Input SMILES string
compound_name = input("Name: ")
smiles = input("Smiles: ")
smiles = remove_stars(smiles)
smiles_fix = duplicate_ring_number_handler(smiles)
mol = Chem.MolFromSmiles(smiles)
atom_properties = get_atom_properties(mol)
prop_copy = copy.deepcopy(atom_properties)
final = ["" for _ in atom_properties]
conn_mat = connectivity_matrix(mol, len(atom_properties))
bead_dict = get_m3_dict()
smiles_token = parse_smiles(smiles_fix)
mapping = map_molecule(smiles_token, conn_mat, atom_properties)
mapping_copy = copy.deepcopy(mapping)
print(mapping_copy)
final = map_martini_beads(mapping, final, bead_dict)
beads_data, connections_data = export_bead_mapping(final, mapping_copy, smiles, compound_name)
gro_maker(beads_data, compound_name)
itp_maker(beads_data, connections_data, compound_name)