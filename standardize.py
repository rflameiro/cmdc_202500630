import numpy as np
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from chembl_structure_pipeline import standardizer


def new_standardize_smiles(smiles, return_mol=True):
    """
    Input = SMILES
    return_mol: Return RDKit Mol object (True) or SMILES (False)
    """
    try:
        molecule = Chem.MolFromSmiles(smiles)
    except:
        return
    
    try:
        # Apply standardize_mol
        std1_mol = standardizer.standardize_mol(molecule)
        # Select largest fragment
        desalter = rdMolStandardize.LargestFragmentChooser()
        desalt_mol = desalter.choose(std1_mol)
        # Apply standardize_mol to the largest fragment
        std2_mol = standardizer.standardize_mol(desalt_mol)
        # Select tautomer with highest score
        te = rdMolStandardize.TautomerEnumerator()
        curated_mol = te.Canonicalize(std2_mol)
        
        if not curated_mol:
            return
        
        if return_mol:
            return curated_mol
        else:
            return Chem.MolToSmiles(curated_mol)
    
    except Exception as e:
        return

# Check for compounds with non-organic atoms
not_organic_pat = Chem.MolFromSmarts("[!#5;!#6;!#7;!#8;!#16;!#15;!F;!Cl;!Br;!I;!#1]")

def non_organic(smiles):
    # Returns True if a non-organic atom is found 
    # or SMILES failed conversion to RDKit Mol
    try:
        mol = Chem.MolFromSmiles(smiles)
        return bool(mol.GetSubstructMatch(not_organic_pat))
    except:  # failed to convert to rdkit_mol
        return True