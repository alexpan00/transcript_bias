import pandas as pd
import sys
import os

def main():
    if len(sys.argv) < 3:
        print("Usage: python merge_oarfish.py <output_file> <input_file1> <input_file2> ...")
        sys.exit(1)
        
    # output csv
    output_file = sys.argv[1]
    
    # read the input files
    input_files = sys.argv[2:]
    
    # read the first df
    sample_name = os.path.basename(input_files[0]).split('.')[0]
    df = pd.read_csv(input_files[0], sep='\t')
    id_col = 'tname' if 'tname' in df.columns else df.columns[0]
    num_col = 'num_reads' if 'num_reads' in df.columns else df.columns[1]
    
    counts_df = df[[id_col, num_col]].rename(columns={num_col: sample_name, id_col: 'tname'})
    
    # loop over the other files and outer merge on tname
    for input_file in input_files[1:]:
        sample_name = os.path.basename(input_file).split('.')[0]
        df2 = pd.read_csv(input_file, sep='\t')
        id_col_2 = 'tname' if 'tname' in df2.columns else df2.columns[0]
        num_col_2 = 'num_reads' if 'num_reads' in df2.columns else df2.columns[1]
        df2_counts = df2[[id_col_2, num_col_2]].rename(columns={num_col_2: sample_name, id_col_2: 'tname'})
        counts_df = pd.merge(counts_df, df2_counts, on='tname', how='outer')
    
    counts_df.fillna(0, inplace=True)
    # write the output file
    counts_df.to_csv(output_file, sep='\t', index=False)
    
if __name__ == '__main__':
    main()