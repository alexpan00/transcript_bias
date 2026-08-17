import argparse
import os

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
    current_transcript_id = None
    current_dedup_id = None
    
    with open(args.input_gtf, "r") as infile, open(args.output_gtf, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                continue
            
            fields = line.strip().split("\t")
            if len(fields) < 9:
                outfile.write(line)
                continue

            feature_type = fields[2]
            attributes = fields[8]
            
            # Simple extractor for transcript_id attribute
            # Looks for: transcript_id "ID";
            import re
            match = re.search(r'transcript_id\s+"([^"]+)"', attributes)
            transcript_id = match.group(1) if match else None

            if feature_type == "transcript":
                current_transcript_id = transcript_id
                
                # Check for duplicate
                new_id = transcript_id
                if new_id in seen_names:
                    i = 1
                    while f"{new_id}_{i}" in seen_names:
                        i += 1
                    new_id = f"{new_id}_{i}"
                
                seen_names.add(new_id)
                current_dedup_id = new_id
                
                # Update attribute string
                if transcript_id:
                    new_attributes = attributes.replace(f'transcript_id "{transcript_id}"', f'transcript_id "{new_id}"')
                    fields[8] = new_attributes
            
            elif feature_type == "exon":
                # Propagate deduplicated ID if this exon belongs to the current transcript
                if transcript_id and transcript_id == current_transcript_id and current_dedup_id:
                    new_attributes = attributes.replace(f'transcript_id "{transcript_id}"', f'transcript_id "{current_dedup_id}"')
                    fields[8] = new_attributes

            outfile.write("\t".join(fields) + "\n")

    print(f"Deduplicated GTF written to {args.output_gtf}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Deduplicate transcript_id in a GTF file")
    add_args(parser)
    args = parser.parse_args()
    main(args)