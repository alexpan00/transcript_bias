#!/usr/bin/env python3
from scipy.io import mmread
import pandas as pd
import argparse

def main(matrix, metadata, transcripts):
    COOmatrix = mmread(matrix)
    # transpose the matrix
    adata = pd.DataFrame.sparse.from_spmatrix(COOmatrix).T
    
    cols = pd.read_csv(metadata)["sample"]
    rows = pd.read_csv(transcripts, header=None, names=["transcript_id"], sep="\t")
    
    # set the columns and rows names of adata
    adata.columns = cols
    adata.index = rows["transcript_id"]

    return adata

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Converts kallisto quant to tsv")
    p.add_argument('-c', action="store", dest="counts", help="matrix abundance file" )
    p.add_argument('-m', action="store", dest="metadata", help="metadata file" )
    p.add_argument('-t', action="store", dest="transcripts", help="transcripts file" )

    p.add_argument("--output", help="path to save tsv counts", default="kallisto_counts.tsv")

    args = p.parse_args()

    counts = main(args.counts, args.metadata, args.transcripts)
    counts.to_csv(args.output, sep="\t")