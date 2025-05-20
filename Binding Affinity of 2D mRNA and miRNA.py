import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import RNA
import pandas as pd
import uuid
import argparse
import sys
import os

class MiRNATargetPredictor:
    def __init__(self, miRNA_seq, mRNA_seq, miRNA_id, mRNA_id):
        self.miRNA_seq = miRNA_seq.upper().replace('T', 'U')
        self.mRNA_seq = mRNA_seq.upper().replace('T', 'U')
        self.miRNA_id = miRNA_id
        self.mRNA_id = mRNA_id
        self.seed_types = {
            '8mer': (2, 8, True),
            '7mer_m8': (2, 8, False),
            '7mer_A1': (2, 7, True),
            '6mer_2-7': (2, 7, False),
            '6mer_3-8': (3, 8, False)
        }
        
    def reverse_complement(self, seq):
        """Generate reverse complement of RNA sequence"""
        return str(Seq(seq).complement_rna())[::-1]
    
    def check_seed_match(self, miRNA_seed, mRNA_segment, seed_type):
        """Check for seed matching between miRNA and mRNA"""
        start, end, needs_A = self.seed_types[seed_type]
        miRNA_segment = miRNA_seed[start-1:end]
        mRNA_rc = self.reverse_complement(mRNA_segment)
        
        # Check Watson-Crick pairing
        pairs = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
        match = all(pairs.get(miRNA_segment[i], '') == mRNA_rc[i] 
                   for i in range(len(miRNA_segment)))
        
        # Check for A at position 1 if required
        if needs_A:
            return match and mRNA_segment[0] == 'A'
        return match
    
    def calculate_free_energy(self, miRNA_seq, mRNA_seq):
        """Calculate free energy of miRNA-mRNA duplex using RNAcofold"""
        try:
            duplex = f"{miRNA_seq}&{mRNA_seq}"
            fc = RNA.fold_compound(duplex)
            _, mfe = fc.mfe()
            return mfe
        except Exception:
            return float('inf')  # Return invalid energy if calculation fails
    
    def assess_site_accessibility(self, mRNA_seq, start, end):
        """Assess structural accessibility of target site"""
        try:
            fc = RNA.fold_compound(mRNA_seq)
            _, mfe = fc.mfe()
            
            # Calculate accessibility by comparing folded and unfolded energies
            subseq = mRNA_seq[max(0, start-10):end+10]
            fc_sub = RNA.fold_compound(subseq)
            _, mfe_sub = fc_sub.mfe()
            
            return mfe - mfe_sub
        except Exception:
            return float('inf')  # Return invalid accessibility if calculation fails
    
    def find_target_sites(self, energy_threshold=-10, window_size=30):
        """Find potential miRNA target sites in mRNA"""
        results = []
        
        for i in range(len(self.mRNA_seq) - window_size + 1):
            mRNA_window = self.mRNA_seq[i:i+window_size]
            
            # Check all seed types
            for seed_type in self.seed_types:
                if len(mRNA_window) >= 8:  # Ensure window is long enough
                    if self.check_seed_match(self.miRNA_seq, mRNA_window, seed_type):
                        # Calculate free energy
                        energy = self.calculate_free_energy(self.miRNA_seq, mRNA_window)
                        
                        if energy <= energy_threshold:
                            # Assess site accessibility
                            accessibility = self.assess_site_accessibility(
                                self.mRNA_seq, i, i+window_size)
                            
                            results.append({
                                'miRNA_id': self.miRNA_id,
                                'mRNA_id': self.mRNA_id,
                                'position': i,
                                'seed_type': seed_type,
                                'free_energy': energy,
                                'accessibility': accessibility,
                                'mRNA_sequence': mRNA_window
                            })
        
        return pd.DataFrame(results)
    
    def rank_targets(self, df):
        """Rank potential targets based on combined scores"""
        if df.empty:
            return df
            
        # Normalize scores
        df['energy_score'] = (df['free_energy'] - df['free_energy'].min()) / \
                           (df['free_energy'].max() - df['free_energy'].min() + 1e-10)
        df['access_score'] = (df['accessibility'] - df['accessibility'].min()) / \
                           (df['accessibility'].max() - df['accessibility'].min() + 1e-10)
        
        # Combine scores (weighted average)
        df['combined_score'] = 0.6 * df['energy_score'] + 0.4 * df['access_score']
        
        # Sort by combined score
        return df.sort_values('combined_score', ascending=False)

def get_valid_text_file_path(prompt):
    """Prompt user for a valid .txt file path until one is provided or user exits"""
    while True:
        file_path = input(prompt).strip()
        if file_path.lower() in ['exit', 'quit']:
            print("Exiting program.")
            sys.exit(0)
        if not file_path.lower().endswith('.txt'):
            print(f"Error: File '{file_path}' must be a .txt file. Please try again or type 'exit' to quit.")
            continue
        if os.path.isfile(file_path) and os.access(file_path, os.R_OK):
            return file_path
        else:
            print(f"Error: File '{file_path}' does not exist or is not readable. Please try again or type 'exit' to quit.")

def validate_text_file(input_file):
    """Validate that the input text file is tab-separated and has required columns"""
    try:
        # Read the first few lines to check format
        with open(input_file, 'r') as f:
            first_line = f.readline().strip()
            if not first_line:
                raise ValueError("Input file is empty.")
            if '\t' not in first_line:
                raise ValueError("Input file must be tab-separated.")
        
        # Read with pandas to check columns
        df_input = pd.read_csv(input_file, sep='\t', nrows=1)
        required_columns = ['mirna_id', 'mirna_seq', 'mrna_id', 'mrna_seq', 'label']
        if not all(col in df_input.columns for col in required_columns):
            raise ValueError(
                f"Input file must contain columns: {', '.join(required_columns)}. "
                f"Found: {', '.join(df_input.columns)}"
            )
        return True
    except Exception as e:
        raise ValueError(f"Invalid input file format: {str(e)}")

def process_input_file(input_file, output_file):
    """Process input text file and predict miRNA targets"""
    try:
        # Validate text file format
        validate_text_file(input_file)
        
        # Read input file
        df_input = pd.read_csv(input_file, sep='\t')
        
        all_results = []
        
        # Process each miRNA-mRNA pair
        for _, row in df_input.iterrows():
            predictor = MiRNATargetPredictor(
                miRNA_seq=row['mirna_seq'],
                mRNA_seq=row['mrna_seq'],
                miRNA_id=row['mirna_id'],
                mRNA_id=row['mrna_id']
            )
            targets = predictor.find_target_sites()
            ranked_targets = predictor.rank_targets(targets)
            
            if not ranked_targets.empty:
                all_results.append(ranked_targets)
        
        # Combine all results
        if all_results:
            final_df = pd.concat(all_results, ignore_index=True)
            # Select relevant columns for output
            output_columns = ['miRNA_id', 'mRNA_id', 'position', 'seed_type', 
                           'free_energy', 'accessibility', 'combined_score', 'mRNA_sequence']
            final_df = final_df[output_columns]
            
            # Save to Excel
            final_df.to_excel(output_file, index=False, engine='openpyxl')
            print(f"Results saved to {output_file}")
        else:
            print("No valid target sites found for any miRNA-mRNA pair.")
            # Create empty Excel file
            pd.DataFrame(columns=['miRNA_id', 'mRNA_id', 'position', 'seed_type', 
                                'free_energy', 'accessibility', 'combined_score', 'mRNA_sequence'])\
                .to_excel(output_file, index=False, engine='openpyxl')
            print(f"Empty results saved to {output_file}")
            
    except Exception as e:
        print(f"Error processing input file: {str(e)}")
        sys.exit(1)

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Predict miRNA target sites from a tab-separated text file")
    parser.add_argument('--input', help="Input text file path (.txt, tab-separated)")
    parser.add_argument('--output', default='miRNA_targets.xlsx', 
                       help="Output Excel file path (default: miRNA_targets.xlsx)")
    
    args = parser.parse_args()
    
    # Get input file path
    input_file = args.input
    if input_file:
        if not input_file.lower().endswith('.txt'):
            print(f"Error: Input file '{input_file}' must be a .txt file.")
            sys.exit(1)
        if not (os.path.isfile(input_file) and os.access(input_file, os.R_OK)):
            print(f"Error: Input file '{input_file}' does not exist or is not readable.")
            sys.exit(1)
        try:
            validate_text_file(input_file)
        except ValueError as e:
            print(f"Error: {str(e)}")
            sys.exit(1)
    else:
        input_file = get_valid_text_file_path(
            "Please enter the path to the input text file (.txt, tab-separated, with columns: mirna_id, mirna_seq, mrna_id, mrna_seq, label)\n"
            "or type 'exit' to quit: "
        )
        try:
            validate_text_file(input_file)
        except ValueError as e:
            print(f"Error: {str(e)}")
            sys.exit(1)
    
    # Process the input file
    process_input_file(input_file, args.output)

if __name__ == "__main__":
    main()