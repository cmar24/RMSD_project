import os
import pandas as pd
from Bio.PDB import PDBParser, Superimposer
import numpy as np
import argparse

class CustomCGRMSD:
    def __init__(self, atom_names):
        self.atom_names = atom_names
        self.parser = PDBParser(QUIET=True)
        self.rmsd_column_name = f"CGRMSD_{'_'.join(atom_names)}"
    
    def get_atoms(self, structure):
        """Retrieve specific atoms from a structure."""
        atoms = []
        for residue in structure.get_residues():
            for atom in residue:
                if atom.get_id() in self.atom_names:
                    atoms.append(atom)
        return atoms
    
    def calculate_rmsd(self, native_atoms, predict_atoms):
        """Calculate RMSD for aligned atoms."""
        diffs = np.array([native.coord - predict.coord for native, predict in zip(native_atoms, predict_atoms)])
        rmsd = np.sqrt((diffs ** 2).sum() / len(diffs))
        return rmsd
    
    def process_single_structure(self, native_file, predict_file):
        """Calculate CGRMSD for a single native and predicted structure."""
        # Load native structure
        native_structure = self.parser.get_structure('native', native_file)
        native_atoms = self.get_atoms(native_structure)
        
        # Load predicted structure
        predict_structure = self.parser.get_structure('predict', predict_file)
        predict_atoms = self.get_atoms(predict_structure)
        
        # Check if atom counts match
        if len(native_atoms) != len(predict_atoms):
            raise ValueError("The number of atoms in the native and predicted structures do not match.")
        
        # Superimpose and calculate RMSD
        super_imposer = Superimposer()
        super_imposer.set_atoms(native_atoms, predict_atoms)
        super_imposer.apply(predict_structure.get_atoms())
        
        rmsd = self.calculate_rmsd(native_atoms, predict_atoms)
        return rmsd

    def process_all_structures(self, native_dir, predict_dir, csv_dir):
        """Process all structures in the native directory."""
        for native_file in os.listdir(native_dir):
            structure_name = native_file.split('.')[0]
            self.process_structure(structure_name, native_dir, predict_dir, csv_dir)

    def process_structure(self, structure_name, native_dir, predict_dir, csv_dir):
        """Calculate and save CGRMSD for each predicted structure of a given structure."""
        # Load native structure
        native_path = os.path.join(native_dir, f"{structure_name}.pdb")
        native_structure = self.parser.get_structure(structure_name, native_path)
        native_atoms = self.get_atoms(native_structure)
        
        # Load CSV file for scores and predicted structures
        csv_path = os.path.join(csv_dir, f"{structure_name}.csv")
        df = pd.read_csv(csv_path, dtype={0: str})  # Ensure first column is treated as string
        rmsd_values = [] 
        
        # Calculate CGRMSD for each predicted structure
        for _, row in df.iterrows():
            predict_file = str(row.iloc[0]).replace("normalized_", "")  # Ensure predict_file is string
            predict_path = os.path.join(predict_dir, structure_name, predict_file)

            if not os.path.exists(predict_path):
                rmsd_values.append(None)
                continue
            
            # Load predicted structure
            predict_structure = self.parser.get_structure(f"{structure_name}_pred", predict_path)
            predict_atoms = self.get_atoms(predict_structure)
            
            # Check if atom counts match
            if len(native_atoms) != len(predict_atoms):
                rmsd_values.append(None)
                continue
            
            # Superimpose and calculate RMSD
            super_imposer = Superimposer()
            super_imposer.set_atoms(native_atoms, predict_atoms)
            super_imposer.apply(predict_structure.get_atoms())
            
            rmsd = self.calculate_rmsd(native_atoms, predict_atoms)
            rmsd_values.append(rmsd)
        
        # Save updated CSV with the custom column name for CGRMSD
        df[self.rmsd_column_name] = rmsd_values
        df.to_csv(csv_path, index=False)

    def calculate_correlation(self, df):
        """Calculate and print the Pearson correlation matrix for a given DataFrame."""
        # Drop rows with missing values
        df_filtered = df.dropna()     

        # Select only numeric columns
        df_numeric = df_filtered.select_dtypes(include=[np.number])

        # Calculate correlation matrix using Pearson method
        correlation_matrix = df_numeric.corr(method='pearson')

        return correlation_matrix

def main():
    parser = argparse.ArgumentParser(description="Calculate CGRMSD and correlation matrix for RNA structures.")
    parser.add_argument('input_type', type=str, choices=['single', 'all'], help="Specify whether to process 'single' pair of files or 'all' structures in directories.")
    parser.add_argument('atom_names', type=str, nargs='+', help="List of atom names to be considered for CGRMSD calculation")
    
    # Arguments for 'single' input_type
    parser.add_argument('--native_file', type=str, help="Path to the native PDB file")
    parser.add_argument('--predict_file', type=str, help="Path to the predicted PDB file")
    parser.add_argument('--csv_file', type=str, help="Path to the CSV file containing scores and predicted structures")
    
    # Arguments for 'all' input_type
    parser.add_argument('--native_dir', type=str, help="Path to the directory containing native PDB files")
    parser.add_argument('--predict_dir', type=str, help="Path to the directory containing predicted PDB files")
    parser.add_argument('--csv_dir', type=str, help="Path to the directory containing CSV files")
    
    args = parser.parse_args()
    
    cgrmsd_calculator = CustomCGRMSD(args.atom_names)
    
    if args.input_type == 'single':
        if not args.native_file or not args.predict_file or not args.csv_file:
            print("For 'single' input_type, --native_file, --predict_file, and --csv_file must be specified.")
            return
        
        # Calculate and print the CGRMSD
        rmsd = cgrmsd_calculator.process_single_structure(args.native_file, args.predict_file)
        print(f"The calculated CGRMSD is: {rmsd}")
        
        # Calculate and print the Pearson correlation matrix
        # Ommited because it does not make sense to calculate the correlation matrix here due to the lack of variation
        #df = pd.read_csv(args.csv_file)
        #predict_file_basename = os.path.basename(args.predict_file)
        #row = df[df.iloc[:, 0].str.replace("normalized_", "") == predict_file_basename]
        # Check if the row exists
        #if row.empty:
        ##    print(f"No matching entry found for {predict_file_basename} in the CSV file.")
        #    return
        # Add a new column for the calculated CGRMSD in the extracted row
        #row.loc[:, cgrmsd_calculator.rmsd_column_name] = rmsd
        #print(row)
        #correlation_matrix = cgrmsd_calculator.calculate_correlation(row)
        # Print the correlation matrix
        #print("Pearson correlation matrix:")
        #print(correlation_matrix)
    
    elif args.input_type == 'all':
        if not args.native_dir or not args.predict_dir or not args.csv_dir:
            print("For 'all' input_type, --native_dir, --predict_dir, and --csv_dir must be specified.")
            return
        
        # Process all structures and calculate CGRMSD
        cgrmsd_calculator.process_all_structures(args.native_dir, args.predict_dir, args.csv_dir)
        correlation_matrix_dir = os.path.join(args.csv_dir, "cormat")
        os.makedirs(correlation_matrix_dir, exist_ok=True)

        # Calculate and save the Pearson correlation matrix for each structure
        for csv_file in os.listdir(args.csv_dir):
            if not csv_file.endswith('.csv'):
                continue  # Skip non-CSV files
            protein_name = csv_file.split('.')[0]
            csv_path = os.path.join(args.csv_dir, csv_file)
            df = pd.read_csv(csv_path)
            correlation_matrix = cgrmsd_calculator.calculate_correlation(df)
            correlation_matrix_path = os.path.join(correlation_matrix_dir, f"{protein_name}_cormat.csv")
            correlation_matrix.to_csv(correlation_matrix_path)


if __name__ == "__main__":
    main()

