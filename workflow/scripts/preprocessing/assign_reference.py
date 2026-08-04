# Assign isoforms to the reference transcripts when they are reference match.
# In case that multiple isoforms match the same reference transcript sum the
# expression of the isoforms.


import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="""
    Collapse transcript quantifications based on sqanti classification.

    Transcripts with subcategory == 'reference_match' and the same associated_transcript
    will be collapsed (counts summed), and the first transcript's info will be retained.
    """)
    parser.add_argument("-c", "--classification", required=True, help="Path to SQANTI3 classification file")
    parser.add_argument("-q", "--quantifications", required=True, help="Path to quantification matrix CSV (transcript_id as index)")
    parser.add_argument("-m", "--metadata", required=True, help="Path to transcript stats file")
    parser.add_argument("-o", "--output_prefix", default="collapsed", help="Prefix for output files")
    return parser.parse_args()


def collapse_transcripts(meta:pd.DataFrame, quant:pd.DataFrame, stats: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame,pd.DataFrame]:
    '''
    Given a classification file, a quantification matrix, and a stats file,
    assign isoforms to the reference transcripts when they are reference match.
    In case that multiple isoforms match the same reference transcript sum the
    expression of the isoforms and compute the median length and GC content.
    The first isoform's info will be retained.

    Args:
        meta (pd.DataFrame): SQANTI3 classification file
        quant (pd.DataFrame): quantification matrix
        stats (pd.DataFrame): transcript id and GC content

    Returns:
        pd.Dataframe: returns 3 dataframes similar to the inputs but with the
        information collapsed by reference transcript.
        The 4th dataframe is a mapping of the collapsed transcripts to the
        original transcripts.
    '''
    # Filter transcripts to collapse
    # Adjust subcategory for mono-exons that match reference criteria
    mono_mask = (
        (meta['exons'] == 1) & 
        (meta['structural_category'] == 'full-splice_match') & 
        (meta['ref_exons'] == 1) & 
        (meta['diff_to_TSS'].abs() < 50) & 
        (meta['diff_to_TTS'].abs() < 50)
    )
    meta.loc[mono_mask, 'subcategory'] = 'reference_match'

    filtered = meta[meta['subcategory'] == 'reference_match'].copy()
    grouped = filtered.groupby('associated_transcript')

    new_meta_rows = []
    new_quant_rows = []
    new_stats_rows = []
    collapse_map_rows = []
    # for all the isoforms that match the same reference transcript as FSM
    # reference-match, sum their counts and take the median of their GC and length.
    # and keep the info of the classification of the first isoform in the group.
    # Replace the isoform id with the associated transcript id (reference transcript id)
    for assoc_tr, group in grouped:
        transcripts = group['isoform'].tolist()
        valid_quant_tr = [t for t in transcripts if t in quant.index]
        valid_stats_tr = [t for t in transcripts if t in stats.index]

        summed_counts = quant.loc[valid_quant_tr].sum() if valid_quant_tr else pd.Series(0, index=quant.columns)
        grouped_stats = stats.loc[valid_stats_tr].median(numeric_only=True) if valid_stats_tr else pd.Series(dtype=float)

        # Metadata from first entry
        rep_row = group.iloc[0].copy()
        rep_row['isoform'] = assoc_tr
        new_meta_rows.append(rep_row)

        # Summed quant row
        summed_counts.name = assoc_tr
        new_quant_rows.append(summed_counts)
        
        # Stats from median computation
        grouped_stats.name = assoc_tr
        new_stats_rows.append(grouped_stats)

        # Collapse map
        collapse_map_rows.append({
            'new_transcript_id': assoc_tr,
            'collapsed_transcripts': ",".join(transcripts)
        })

    # Create result DataFrames with explicit column structures
    new_meta = pd.DataFrame(new_meta_rows, columns=meta.columns) if new_meta_rows else pd.DataFrame(columns=meta.columns)
    new_quant = pd.DataFrame(new_quant_rows, columns=quant.columns) if new_quant_rows else pd.DataFrame(columns=quant.columns)
    new_stats = pd.DataFrame(new_stats_rows, columns=stats.columns) if new_stats_rows else pd.DataFrame(columns=stats.columns)
    collapse_map = pd.DataFrame(collapse_map_rows, columns=['new_transcript_id', 'collapsed_transcripts']) if collapse_map_rows else pd.DataFrame(columns=['new_transcript_id', 'collapsed_transcripts'])

    # Add back non-collapsed transcripts
    collapsed_ids = set(filtered['isoform'])
    remaining_meta = meta[~meta['isoform'].isin(collapsed_ids)].copy()
    valid_rem_quant = [t for t in remaining_meta['isoform'] if t in quant.index]
    valid_rem_stats = [t for t in remaining_meta['isoform'] if t in stats.index]
    remaining_quant = quant.loc[valid_rem_quant]
    remaining_stats = stats.loc[valid_rem_stats]

    # Combine everything
    final_meta = pd.concat([new_meta, remaining_meta], ignore_index=True)
    final_quant = pd.concat([new_quant, remaining_quant])
    final_stats = pd.concat([new_stats, remaining_stats])

    return final_meta, final_quant, final_stats, collapse_map


def main():
    args = parse_args()

    # Load input
    meta = pd.read_csv(args.classification, sep="\t")
    quant = pd.read_csv(args.quantifications, sep="\t", index_col=0)
    stats = pd.read_csv(args.metadata, sep="\t", index_col=0)

    # Run collapsing
    final_meta, final_quant, final_stats, collapse_map = collapse_transcripts(meta, quant, stats)

    # Save outputs
    final_meta.to_csv(f"{args.output_prefix}_classification.txt", index=False, sep="\t")
    final_quant.to_csv(f"{args.output_prefix}_counts.tsv", sep="\t", index_label="superPBID")
    final_stats.to_csv(f"{args.output_prefix}_transcript_models.tsv", sep="\t", index_label="#id")
    collapse_map.to_csv(f"{args.output_prefix}_collapse_map.csv", index=False)
    
    # Create final transcript 2 gene map based on final_meta (collapsed sqanti classification)
    tx2gene = final_meta[['isoform', 'associated_gene']].copy()
    # Don't write header and index
    tx2gene.to_csv(f"{args.output_prefix}_transcript_to_gene.tsv", sep="\t", index=False, header=False)

if __name__ == "__main__":
    main()
