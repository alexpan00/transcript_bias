import argparse
from CAGE_class import CAGEPeak
from GTF_classes import read_gtf_as_transcripts, Transcript, GTFExon
from collections import defaultdict
from joblib import Parallel, delayed
import pandas as pd
from os.path import join
from os import makedirs
import os


def get_matches(trans_dict, cage_dict):
    '''
    Function to get candidate transcripts for collapsing based on CAGE data.
    '''
    cage_hits = list()
    no_hits = set()
    for trans_id, tx in trans_dict.items():
        start = tx.start - 1 if tx.strand == '+' else tx.end -1 # BED is 0-based and GTF is 1-based
        hit = cage_dict[(tx.seqname, tx.strand)].find(start, start + 1)
        if hit:
            cage_hits.append([tx.attributes["transcript_id"], hit[0], tx.get_UJC(), tx.start, tx.end])
        else:
            no_hits.add(trans_id)
    return (cage_hits, no_hits)


def cluster_transcripts(common_UJC: list, start: int, end: int):
    '''
    Given a list of transcripts that have the same UJC, cluster them based on their
    junctions. We just select the longest end for the UJC.

    Args:
        common_UJC (list): list of transcripts that have the same UJC.
        
    Returns:
        transcript (Transcript): A new transcript that has the same UJC and the 
                                    longest ends.
    '''
    
    # Get the attributes of the first transcript
    strand = common_UJC[0].strand
    chrom = common_UJC[0].seqname
    source = common_UJC[0].source
    score = common_UJC[0].score
    attr_str = common_UJC[0].attr_str()

    # Create a new transcript with the same UJC
    transcript = Transcript(chrom, source, "transcript", start, end, score, strand, ".", attr_str)
    num_exons = len(common_UJC[0].exons)
    for exon_idx, exon in enumerate(common_UJC[0].exons):
        e_start = start if exon_idx == 0 else exon.start
        e_end = end if exon_idx == num_exons - 1 else exon.end
        new_exon = GTFExon(chrom, source, "exon", e_start, e_end, score, strand, ".", attr_str)
        transcript.add_exon(new_exon)
    return transcript
        

def collapse_matches(hits: list, start: int, end: int, transcript_dict):
    '''
    Given a list of transcript that have a their TSS in the same CAGE peak, collapse
    the ones that have the same UJC 

    Args:
        hits (list): A lsit of transcripts that have the same CAGE peak and UJC.
        start (int): The start of the collapsed transcript
        end (int): The end of the collapsed transcript
        transcript_dict (dict): A dictionary of transcripts.
    '''
    
    # Collapse the transcripts
    if len(hits) > 1:
        # Get the longest transcript
        final_transcript = cluster_transcripts([transcript_dict[hit] for hit in hits], start, end)
        trans_cluster_info = {id1: final_transcript.attributes["transcript_id"] for id1 in hits[1:]}
    else:
        final_transcript = transcript_dict[hits[0]]
        trans_cluster_info = {}
    return final_transcript, trans_cluster_info
            
def group_transcripts_by_chr_strand(trans_dict):
    grouped = defaultdict(dict)
    for tid, tx in trans_dict.items():
        key = (tx.seqname, tx.strand)
        grouped[key][tid] = tx
    return grouped

def main():
    parser = argparse.ArgumentParser(description="Collapse transcript models based on CAGE data.")
    parser.add_argument("cage_file", help="Input file containing CAGE data.")
    parser.add_argument("gtf_file", help="GTF file containing transcript models.")
    parser.add_argument("counts_file", help="Coutns file.")
    parser.add_argument("output_prefix", help="Output directory for collapsed transcripts.")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use for parallel processing.")
    
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    out_dir, out_name = os.path.split(args.output_prefix)
    if out_dir:
        makedirs(out_dir, exist_ok=True)
    
    
    print(f"Collapsing {args.gtf_file} based on CAGE data from {args.cage_file}")
    # get final final files names
    out_gtf = join(out_dir, f"{out_name}.cage.gtf")
    out_counts = join(out_dir, f"{out_name}.cage_counts.txt")
    out_report = join(out_dir, f"{out_name}.cage_report.txt")
    
    
    # Read CAGE data    
    cage = CAGEPeak(args.cage_file)
    
    # Read counts file
    counts = pd.read_csv(args.counts_file, index_col=0)

    # Read GTF file
    trans_dict = read_gtf_as_transcripts(args.gtf_file)
    # Group transcripts by chromosome and strand
    grouped_transcripts = group_transcripts_by_chr_strand(trans_dict)

    # Get CAGE hits for each transcript
    results = Parallel(n_jobs=args.threads, backend='threading')(delayed(get_matches)(grouped_transcripts[key], cage.cage_peaks) for key in grouped_transcripts.keys())
    all_cage_res = []
    all_no_hits = set()
    for cage_hits, no_hits in results:
        all_cage_res.extend(cage_hits)
        all_no_hits.update(no_hits)
    no_hits_transcripts = [trans_dict[transcript] for transcript in all_no_hits] # get transcript object from ids
    
    cage_hit_df = pd.DataFrame(all_cage_res, columns=["transcript_id", "cage_peak", "UJC", "start", "end"])
    
    # Get the start and end of the transcripts that share CAGE hit and UJC
    res = cage_hit_df.groupby(["UJC", "cage_peak"]).agg({"transcript_id": lambda x: list(x), "start": "min", "end": "max"}).reset_index()
    
    # Generate a merged transcript and cluster information
    results = Parallel(n_jobs=args.threads, backend='threading')(delayed(collapse_matches)(row["transcript_id"], row["start"], row["end"], trans_dict) for _, row in res.iterrows())

    hits_transcripts = []
    with open(out_report, "w") as f:
        f.write("query,target\n")
        cluster_updates = []
        for transcript, cluster in results:
            if cluster:
                cluster_updates.extend(cluster.items())
            hits_transcripts.append(transcript)

        # Batch update counts and write cluster information
        for query, target in cluster_updates:
            f.write(f"{query},{target}\n")
        
        # Aggregate counts for all clusters safely
        for query, target in cluster_updates:
            if query in counts.index and target in counts.index:
                counts.loc[target] += counts.loc[query]
        valid_queries = [query for query, _ in cluster_updates if query in counts.index]
        if valid_queries:
            counts.drop(valid_queries, inplace=True)

    # combine the 2 list of transcripts and sort them
    all_transcripts = hits_transcripts + no_hits_transcripts
    all_transcripts.sort(key=lambda t: (t.seqname, t.start, t.end))
    
    # write the clustered gtf
    with open(out_gtf, "w") as final_gtf:
        for trans in all_transcripts:
            final_gtf.write(str(trans))

    # write the counts file
    counts.to_csv(out_counts, sep=",", index_label="id")

if __name__ == "__main__":
    main()