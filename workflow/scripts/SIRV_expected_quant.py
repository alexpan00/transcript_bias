import argparse
import os
import pandas as pd
import json

parser = argparse.ArgumentParser()
parser.add_argument('--reads', help = 'Number of reads aligned to SIRVs')
parser.add_argument('--set', help = 'Set used in the experiment: E0, E1, E2')
parser.add_argument('--out', help = 'Name of an output document')
parser.add_argument('--mix_matrix', help = 'JSON file with SIRV mix matrix')

args = parser.parse_args()

with open(args.mix_matrix, 'r') as f:
    mix_matrix = json.load(f)

def calculate_expected_counts(reads, set_name, mix_matrix, out):
    sirvs = mix_matrix[set_name]
    count = 0
    for key in sirvs:
        exp_dict = sirvs[key]
        exp = sum(exp_dict.values())
        count += exp
    unit = float(int(reads) / count) if count > 0 else 0.0
    for key in sirvs:
        for iso_key in sirvs[key]:
            row = [iso_key, float(sirvs[key][iso_key] * unit)]
            out.append(row)

out = [['SIRV', 'exp_count']]

if args.reads and os.path.exists(args.reads):
    with open(args.reads, 'r') as f:
        reads = int(f.read().strip())
else:
    reads = int(args.reads) if args.reads else 0

calculate_expected_counts(reads, args.set, mix_matrix, out)

out_df = pd.DataFrame(out)
out_df.to_csv(args.out, header = False, index = False, sep = '\t')
