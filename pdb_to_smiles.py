from rdkit import Chem
from rdkit.Chem import AllChem

def pdb_to_smiles(pdb_file):
    """
    Convert a PDB file to SMILES string
    
    Args:
        pdb_file (str): Path to the PDB file
    
    Returns:
        str: SMILES representation of the molecule
    """
    try:
        # Read molecule from PDB file
        mol = Chem.MolFromPDBFile(pdb_file)
        if mol is None:
            return None
        
        # Generate SMILES
        smiles = Chem.MolToSmiles(mol)
        return smiles
    except Exception as e:
        print(f"Error converting PDB to SMILES: {str(e)}")
        return None

def smiles_to_pdb(smiles, output_file):
    """
    Convert a SMILES string to PDB format with improved error handling
    and multiple embedding attempts
    
    Args:
        smiles (str): SMILES string of the molecule
        output_file (str): Path to save the output PDB file
    
    Returns:
        bool: True if conversion successful, False otherwise
    """
    try:
        # Create molecule from SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print("Failed to create molecule from SMILES")
            return False
        
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Try different embedding methods
        success = False
        
        # Method 1: Standard embedding
        try:
            success = AllChem.EmbedMolecule(mol, randomSeed=42) == 0
        except:
            pass

        # Method 2: Try ETKDG if standard embedding fails
        if not success:
            try:
                parameters = AllChem.ETKDGv3()
                parameters.randomSeed = 42
                success = AllChem.EmbedMolecule(mol, parameters) == 0
            except:
                pass

        # Method 3: Try with increased max iterations
        if not success:
            try:
                parameters = AllChem.ETKDGv3()
                parameters.randomSeed = 42
                parameters.maxIterations = 5000  # Increase max iterations
                success = AllChem.EmbedMolecule(mol, parameters) == 0
            except:
                pass

        if not success:
            print("Failed to generate 3D coordinates")
            return False

        # Try structure optimization
        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except:
            print("Warning: Structure optimization failed, using unoptimized structure")
        
        # Write to PDB file
        writer = Chem.PDBWriter(output_file)
        writer.write(mol)
        writer.close()
        
        return True
    except Exception as e:
        print(f"Error converting SMILES to PDB: {str(e)}")
        return False

def main():
    # Example usage
    pdb_file = "7c2p.EM.monomer11.pdb"
    output_pdb = "output.pdb"
    

    # Convert PDB to SMILES
    print("Testing PDB to SMILES conversion...")
    smiles = pdb_to_smiles(pdb_file)
    if smiles:
        print(f"Generated SMILES: {smiles}")
    else:
        print("Using test SMILES string instead")
    
    # Convert SMILES to PDB
    print("\nTesting SMILES to PDB conversion...")
    if smiles_to_pdb(smiles, output_pdb):
        print(f"Successfully converted SMILES to PDB. Saved as {output_pdb}")
    else:
        print("Failed to convert SMILES to PDB")

if __name__ == "__main__":
    main()