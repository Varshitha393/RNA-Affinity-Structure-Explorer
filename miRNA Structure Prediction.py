import RNA
import os
import csv

def read_input_file(file_path):
    """Read miRNA sequences from a tab-separated file."""
    sequences = []
    try:
        with open(file_path, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            for row in reader:
                if len(row) >= 2:  # Ensure at least name and sequence columns
                    seq_id = row[0].strip()
                    sequence = row[1].strip().upper().replace('T', 'U')
                    if sequence:  # Only process non-empty sequences
                        sequences.append((seq_id, sequence))
        print(f"Successfully read {len(sequences)} miRNA sequences from {file_path}")
    except FileNotFoundError:
        print(f"Error: Input file '{file_path}' not found.")
        return []
    except Exception as e:
        print(f"Error reading file: {str(e)}")
        return []
    return sequences

def predict_structure(sequence, seq_id, rna_type):
    """Predict secondary structure using ViennaRNA's RNAfold and generate SVG."""
    try:
        # Initialize RNAfold model
        md = RNA.md()
        fc = RNA.fold_compound(sequence, md)

        # Predict secondary structure and minimum free energy
        (structure, mfe) = fc.mfe()

        # Generate SVG visualization
        output_file = f"{seq_id}_{rna_type}_structure.svg"
        RNA.svg_rna_plot(sequence, structure, output_file)

        return structure, mfe, output_file
    except Exception as e:
        print(f"Error predicting structure for {seq_id}: {str(e)}")
        return None, None, None

def main(input_file):
    """Main function to process miRNA sequences and predict structures."""
    # Read sequences from input file
    sequences = read_input_file(input_file)
    if not sequences:
        print("No valid miRNA sequences found in the input file.")
        return

    # Process each sequence
    processed_count = 0
    for seq_id, sequence in sequences:
        # Validate sequence
        if not all(base in 'AUCG' for base in sequence):
            print(f"Invalid sequence for {seq_id}: Contains non-AUCG bases.")
            continue

        # Classify as miRNA
        rna_type = "miRNA"

        print(f"\nProcessing {seq_id} ({rna_type})")
        print(f"Sequence: {sequence}")
        print(f"Length: {len(sequence)} nucleotides")

        # Predict structure
        structure, mfe, output_file = predict_structure(sequence, seq_id, rna_type)

        if structure is not None:
            # Output results
            print(f"Secondary Structure: {structure}")
            print(f"Minimum Free Energy: {mfe:.2f} kcal/mol")
            print(f"Structure visualization saved to: {output_file}")
            processed_count += 1
        else:
            print(f"Failed to predict structure for {seq_id}.")

    print(f"\nProcessed {processed_count} of {len(sequences)} miRNA sequences successfully.")

if __name__ == "__main__":
    input_file = "test.txt"
    main(input_file)