import RNA
import os
import csv

def read_input_file(file_path):
    """Read mRNA sequences from a tab-separated file with miRNA and mRNA data."""
    sequences = []
    try:
        with open(file_path, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            header = next(reader)  # Skip header row
            for row in reader:
                if len(row) >= 4:  # Ensure at least mirna_id, mirna_seq, mrna_id, mrna_seq
                    mrna_id = row[2].strip()  # mrna_id is in the third column
                    mrna_seq = row[3].strip().upper().replace('T', 'U')  # mrna_seq is in the fourth column
                    if mrna_seq:  # Only process non-empty sequences
                        sequences.append((mrna_id, mrna_seq))
        print(f"Successfully read {len(sequences)} mRNA sequences from {file_path}")
    except FileNotFoundError:
        print(f"Error: Input file '{file_path}' not found.")
        return []
    except Exception as e:
        print(f"Error reading file: {str(e)}")
        return []
    return sequences

def predict_structure(sequence, seq_id, rna_type):
    """Predict secondary structure for mRNA using ViennaRNA's RNAfold."""
    try:
        # Check sequence length to avoid memory issues
        max_length = 10000  # Reasonable limit for RNAfold
        if len(sequence) > max_length:
            print(f"Sequence {seq_id} is too long ({len(sequence)} nt). Truncating to {max_length} nt.")
            sequence = sequence[:max_length]

        # Initialize RNAfold model with parameters optimized for long mRNA
        md = RNA.md()
        md.max_bp_span = 300  # Increased for longer-range interactions in mRNA
        md.window_size = 500  # Larger window for local folding in long sequences
        fc = RNA.fold_compound(sequence, md)

        # Predict secondary structure and minimum free energy
        (structure, mfe) = fc.mfe()

        # Generate SVG output
        output_file = f"{seq_id}_{rna_type}_structure.svg"
        RNA.svg_rna_plot(sequence, structure, output_file)

        return structure, mfe, output_file
    except Exception as e:
        print(f"Error predicting structure for {seq_id}: {str(e)}")
        return None, None, None

def main(input_file):
    """Main function to process mRNA sequences and predict structures."""
    # Read mRNA sequences from input file
    sequences = read_input_file(input_file)
    if not sequences:
        print("No valid mRNA sequences found in the input file.")
        return

    # Process each sequence
    processed_count = 0
    for seq_id, sequence in sequences:
        # Validate sequence
        if not all(base in 'AUCG' for base in sequence):
            print(f"Invalid sequence for {seq_id}: Contains non-AUCG bases.")
            continue

        # Process as mRNA (no length-based skipping)
        rna_type = "mRNA"
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

    print(f"\nProcessed {processed_count} of {len(sequences)} mRNA sequences successfully.")

if __name__ == "__main__":
    input_file = "test.txt"
    main(input_file)