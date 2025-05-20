import pandas as pd
import sys
import uuid

def read_sequences(file_path):
    """
    Read miRNA and mRNA sequences from a tab-delimited file.
    Expected columns: mirna_id, mirna_seq, mrna_id, mrna_seq, label
    Returns: DataFrame with the data
    """
    try:
        df = pd.read_csv(file_path, sep='\t')
        required_columns = ['mirna_id', 'mirna_seq', 'mrna_id', 'mrna_seq', 'label']
        if not all(col in df.columns for col in required_columns):
            raise ValueError("Input file must contain columns: mirna_id, mirna_seq, mrna_id, mrna_seq, label")

        # Convert sequences to uppercase and replace T with U
        df['mirna_seq'] = df['mirna_seq'].str.upper().str.replace('T', 'U')
        df['mrna_seq'] = df['mrna_seq'].str.upper().str.replace('T', 'U')
        valid_bases = set('AUGC')

        # Validate sequences and collect invalid ones
        invalid_rows = []
        for idx, row in df.iterrows():
            invalid_mirna = [base for base in row['mirna_seq'] if base not in valid_bases]
            invalid_mrna = [base for base in row['mrna_seq'] if base not in valid_bases]
            if invalid_mirna or invalid_mrna:
                invalid_rows.append(
                    f"Index {idx}: miRNA {row['mirna_id']} has invalid bases {invalid_mirna}, "
                    f"mRNA {row['mrna_id']} has invalid bases {invalid_mrna}"
                )

        if invalid_rows:
            raise ValueError("Invalid bases found in sequences:\n" + "\n".join(invalid_rows))

        return df

    except FileNotFoundError:
        raise FileNotFoundError(f"File {file_path} not found.")
    except Exception as e:
        raise Exception(f"Error reading file: {str(e)}")

def calculate_binding_affinity(mirna_seq, mrna_seq):
    """
    Calculate binding affinity score between miRNA and mRNA sequences.
    miRNA_seq: 5' to 3' miRNA sequence
    mrna_seq: 3' to 5' mRNA sequence (reverse complement to align with miRNA)
    Returns: Dictionary with affinity score and alignment details
    """
    # Energy values for base pairs (kcal/mol, simplified)
    energy_pairs = {
        ('A', 'U'): -1.0,
        ('U', 'A'): -1.0,
        ('G', 'C'): -2.0,
        ('C', 'G'): -2.0,
        ('G', 'U'): -0.5,
        ('U', 'G'): -0.5
    }

    # Seed region (positions 2-8 of miRNA, 1-based indexing)
    seed_start, seed_end = 1, 7
    seed_seq = mirna_seq[seed_start:seed_end]

    # Initialize variables
    best_score = float('inf')  # Lower energy is better
    best_alignment = ""
    best_position = -1
    seed_mismatch = float('inf')

    # Slide miRNA over mRNA to find best binding site
    for i in range(len(mrna_seq) - len(mirna_seq) + 1):
        mrna_window = mrna_seq[i:i + len(mirna_seq)]
        score = 0.0
        alignment = []

        # Calculate score for each position
        for j, (mirna_base, mrna_base) in enumerate(zip(mirna_seq, mrna_window)):
            pair = (mirna_base, mrna_base)
            if pair in energy_pairs:
                score += energy_pairs[pair]
                alignment.append(f"{mirna_base}-{mrna_base}")
            else:
                score += 2.0  # Penalty for mismatch
                alignment.append(f"{mirna_base}*{mrna_base}")

        # Check seed region complementarity
        seed_mrna = mrna_window[seed_start:seed_end]
        current_seed_mismatch = sum(1 for m, r in zip(seed_seq, seed_mrna)
                                   if (m, r) not in energy_pairs)
        if current_seed_mismatch > 1:  # Allow max 1 mismatch in seed
            continue

        if score < best_score:
            best_score = score
            best_alignment = alignment
            best_position = i
            seed_mismatch = current_seed_mismatch

    # Convert score to positive scale (higher is better)
    affinity_score = -best_score if best_score != float('inf') else 0.0

    return {
        'affinity_score': affinity_score,
        'binding_position': best_position,
        'alignment': ' '.join(best_alignment) if best_alignment else "No valid alignment",
        'seed_mismatch': seed_mismatch if seed_mismatch != float('inf') else "N/A"
    }

def main():
    try:
        # Hardcode the input and output file paths
        input_file = "miRAW_Test0.txt"  # Specify your input file path here
        output_file = "results.csv"      # Specify your output file path here

        # Read dataset
        df = read_sequences(input_file)

        # Calculate binding affinity for each pair
        results = []
        for idx, row in df.iterrows():
            result = calculate_binding_affinity(row['mirna_seq'], row['mrna_seq'])
            results.append({
                'mirna_id': row['mirna_id'],
                'mrna_id': row['mrna_id'],
                'affinity_score': result['affinity_score'],
                'binding_position': result['binding_position'],
                'alignment': result['alignment'],
                'seed_mismatch': result['seed_mismatch'],
                'label': row['label']
            })

        # Save results to CSV
        result_df = pd.DataFrame(results)
        result_df.to_csv(output_file, index=False)
        print(f"Results saved to {output_file}")

    except Exception as e:
        print(f"Error: {str(e)}")
        return  # Exit gracefully in Jupyter instead of sys.exit

if _name_ == "_main_":
    main()