import argparse
import os
import re

def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--input-gtf",
        required=True,
        help="Input GTF file to process",
    )
    parser.add_argument(
        "--output-gtf",
        default=None,
        help="Output GTF file with deduplicated transcript_id names. Default: input GTF base name with _dedup.gtf suffix",
    )

def main(args: argparse.Namespace) -> None:
    # Prepare output GTF filename
    if args.output_gtf is None:
        base, ext = os.path.splitext(args.input_gtf)
        args.output_gtf = f"{base}_dedup{ext}"

    seen_names = set()
    pattern_tx = re.compile(r'\btranscript_id\s+"([^"]+)"')
    pattern_gene = re.compile(r'\bgene_id\s+"([^"]+)"')
    pattern_exon_num = re.compile(r'\bexon_number\s+"?(\d+)"?')

    current_model_key = None
    current_dedup_id = None
    last_exon_start = None
    last_exon_end = None
    last_exon_number = None

    with open(args.input_gtf, "r") as infile, open(args.output_gtf, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                outfile.write(line)
                continue

            seqname = fields[0]
            feature = fields[2]
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                start = None
                end = None
            strand = fields[6]
            attributes = fields[8]

            match_tx = pattern_tx.search(attributes)
            transcript_id = match_tx.group(1) if match_tx else None

            if transcript_id:
                match_gene = pattern_gene.search(attributes)
                gene_id = match_gene.group(1) if match_gene else None
                model_key = (seqname, strand, gene_id, transcript_id)

                match_exon_num = pattern_exon_num.search(attributes)
                exon_number = int(match_exon_num.group(1)) if match_exon_num else None

                is_explicit_model_line = feature.lower() in ('transcript', 'mrna')
                is_key_change = (model_key != current_model_key)

                is_exon_reset = False
                if not is_explicit_model_line and not is_key_change:
                    if feature.lower() == 'exon':
                        if exon_number is not None and last_exon_number is not None and exon_number <= last_exon_number:
                            is_exon_reset = True
                        elif start is not None and last_exon_start is not None and (start <= last_exon_start or (last_exon_end is not None and start <= last_exon_end)):
                            is_exon_reset = True

                is_new_model = is_explicit_model_line or is_key_change or is_exon_reset

                if is_new_model:
                    new_id = transcript_id
                    if new_id in seen_names:
                        i = 1
                        while f"{transcript_id}_{i}" in seen_names:
                            i += 1
                        new_id = f"{transcript_id}_{i}"
                    seen_names.add(new_id)
                    current_dedup_id = new_id
                    current_model_key = model_key

                    if feature.lower() == 'exon':
                        last_exon_start = start
                        last_exon_end = end
                        last_exon_number = exon_number
                    else:
                        last_exon_start = None
                        last_exon_end = None
                        last_exon_number = None
                else:
                    if feature.lower() == 'exon':
                        last_exon_start = start
                        last_exon_end = end
                        if exon_number is not None:
                            last_exon_number = exon_number

                dedup_id = current_dedup_id
                if dedup_id != transcript_id:
                    fields[8] = re.sub(r'\btranscript_id\s+"' + re.escape(transcript_id) + r'"', f'transcript_id "{dedup_id}"', attributes)

            outfile.write("\t".join(fields) + "\n")

    print(f"Deduplicated GTF written to {args.output_gtf}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Deduplicate transcript_id in a GTF file")
    add_args(parser)
    args = parser.parse_args()
    main(args)