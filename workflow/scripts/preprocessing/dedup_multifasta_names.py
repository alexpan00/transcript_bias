import argparse
import os

def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--input-fasta",
        required=True,
        help="Input multi-FASTA file to process",
    )
    parser.add_argument(
        "--output-fasta",
        default=None,
        help="Output multi-FASTA file with deduplicated names. Default:input fasta base name with _dedup.fasta suffix",
    )

def main(args: argparse.Namespace) -> None:
    # Prepare output FASTA filename
    if args.output_fasta is None:
        base, ext = os.path.splitext(args.input_fasta)
        args.output_fasta = f"{base}_dedup{ext}"

    seen_names = set()
    with open(args.input_fasta, "r") as infile, open(args.output_fasta, "w") as outfile:
        for line in infile:
            if line.startswith(">"):
                orig_name = line[1:].strip()
                name = orig_name
                if orig_name in seen_names:
                    i = 1
                    candidate = f"{orig_name}_{i}"
                    while candidate in seen_names:
                        i += 1
                        candidate = f"{orig_name}_{i}"
                    name = candidate
                    line = f">{name}\n" # rename duplicated transcript name
                seen_names.add(name)
            outfile.write(line)

    print(f"Deduplicated FASTA written to {args.output_fasta}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Deduplicate transcript names in a multi-FASTA file")
    add_args(parser)
    args = parser.parse_args()
    main(args)