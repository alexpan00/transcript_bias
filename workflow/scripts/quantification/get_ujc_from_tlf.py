#!/usr/bin/env python3

import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def parse_attributes(attribute_field: str) -> Dict[str, str]:
    attributes: Dict[str, str] = {}
    for raw_item in attribute_field.strip().split(";"):
        item = raw_item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key.strip()] = value.strip().strip('"').strip("'")
        elif " " in item:
            key, value = item.split(" ", 1)
            attributes[key.strip()] = value.strip().strip('"').strip("'")
    return attributes


def parse_exons(exons_field: str) -> List[Tuple[int, int]]:
    exon_ranges: List[Tuple[int, int]] = []
    for exon_text in exons_field.split(","):
        start_text, end_text = exon_text.split("-", 1)
        exon_ranges.append((int(start_text), int(end_text)))
    return exon_ranges


def build_ujc(chromosome: str, strand: str, exon_ranges: List[Tuple[int, int]]) -> str:
    boundary_parts: List[str] = []
    for exon_index, (exon_start, exon_end) in enumerate(exon_ranges):
        if exon_index > 0:
            boundary_parts.append(str(exon_start))
        if exon_index < len(exon_ranges) - 1:
            boundary_parts.append(str(exon_end))

    return f"{chromosome}_{strand}_" + "_".join(boundary_parts)


def main() -> int:
    if len(sys.argv) != 3:
        print(
            "Usage: get_ujc_from_tlf.py <input.tlf> <output.tsv>",
            file=sys.stderr,
        )
        return 1

    input_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])

    with input_path.open("r", encoding="utf-8") as input_handle, output_path.open(
        "w", encoding="utf-8", newline=""
    ) as output_handle:
        writer = csv.writer(output_handle, delimiter="\t")
        writer.writerow(["transcript_id", "UJC"])
        written = 0
        mono_exonic = 0

        for line_number, raw_line in enumerate(input_handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 9:
                raise ValueError(f"Line {line_number} does not contain 9 TLF fields")

            feature_type = fields[2]
            if feature_type != "transcript":
                continue

            chromosome = fields[0]
            strand = fields[6]
            attributes = parse_attributes(fields[8])
            transcript_id = attributes.get("ID") or attributes.get("transcript_id")
            exons_field = attributes.get("exons")

            if not transcript_id:
                raise ValueError(f"Missing transcript ID on line {line_number}")
            if not exons_field:
                continue

            exon_ranges = parse_exons(exons_field)
            # A UJC is a chain of splice junctions, so a mono-exonic transcript
            # has none to report. They used to share a single
            # "<chrom>_<strand>_mono-exon" identifier, which collapsed every
            # mono-exonic transcript on a chromosome and strand into one UJC.
            if len(exon_ranges) < 2:
                mono_exonic += 1
                continue

            ujc = build_ujc(chromosome, strand, exon_ranges)
            writer.writerow([transcript_id, ujc])
            written += 1

    print(
        f"{written} spliced transcripts written, "
        f"{mono_exonic} mono-exonic transcripts skipped (no splice junctions)"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())