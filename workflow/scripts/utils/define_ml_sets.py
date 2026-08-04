#!/usr/bin/env python3

import pandas as pd
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="Select TP and TN sets from SQANTI3 classification file.")
    parser.add_argument("classification_file", help="Path to SQANTI3 classification file (e.g., *.classification.txt)")
    parser.add_argument("--output_dir", default=".", help="Directory to save TP and TN lists")
    parser.add_argument("--max_size", type=int, default=3000, help="Maximum number of transcripts per set")
    args = parser.parse_args()

    # Load classification file
    df = pd.read_csv(args.classification_file, sep='\t', low_memory=False)

    # Ensure required columns are present
    required_columns = ['structural_category', 'all_canonical', 'within_CAGE_peak', 'within_polyA_site', 'exons', 'isoform']
    for col in required_columns:
        if col not in df.columns:
            raise ValueError(f"Required column '{col}' not found in the classification file.")

    # Helper function for robust boolean parsing
    def parse_bool(s):
        if s.dtype == bool:
            return s
        return s.astype(str).str.upper().isin(['TRUE', '1', 'YES'])

    cage_bool = parse_bool(df['within_CAGE_peak'])
    polya_bool = parse_bool(df['within_polyA_site'])
    exons_num = pd.to_numeric(df['exons'], errors='coerce')

    # Define True Positives (TP)
    tp_conditions = (
        (df['structural_category'] == 'full-splice_match') &
        (df['all_canonical'] == "canonical") &
        (cage_bool == True) &
        (polya_bool == True) &
        (exons_num > 1) 
    )
    tp_df = df[tp_conditions]

    # Define True Negatives (TN)
    tn_conditions = (
        (df['structural_category'] != 'full-splice_match') &
        (
            (df['all_canonical'] == "non_canonical") |
            (cage_bool != True) |
            (polya_bool != True)
        ) &
        (exons_num > 1)
    )
    tn_df = df[tn_conditions]

    # Sample up to max_size transcripts for each set
    tp_sample = tp_df.sample(n=min(len(tp_df), args.max_size), random_state=123)
    tn_sample = tn_df.sample(n=min(len(tn_df), args.max_size), random_state=123)

    # Save isoform IDs to files
    os.makedirs(args.output_dir, exist_ok=True)
    tp_file = os.path.join(args.output_dir, "TP_list.txt")
    tn_file = os.path.join(args.output_dir, "TN_list.txt")
    tp_sample['isoform'].to_csv(tp_file, index=False, header=False)
    tn_sample['isoform'].to_csv(tn_file, index=False, header=False)

    print(f"Saved {len(tp_sample)} TP isoforms to {tp_file}")
    print(f"Saved {len(tn_sample)} TN isoforms to {tn_file}")
    
    # Exclusion file for variables that shouldn't be used in ML filtering
    exclusion_file = os.path.join(args.output_dir, "exclusion_list.txt")
    with open(exclusion_file, 'w') as f:
        f.write("all_canonical\n")
        f.write("within_CAGE_peak\n")
        f.write("within_polyA_site\n")
        f.write("dist_to_CAGE_peak\n")
        f.write("dist_to_polyA_site\n")
    print(f"Exclusion list saved to {exclusion_file}")

if __name__ == "__main__":
    main()
