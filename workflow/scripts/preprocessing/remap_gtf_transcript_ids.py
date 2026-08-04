#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Remap transcript IDs in a GTF/GFF file using a collapse map CSV. "
            "Unmapped transcript IDs are kept unchanged."
        )
    )
    parser.add_argument("input_gtf", help="Input GTF/GFF file path")
    parser.add_argument("map_csv", help="CSV with columns new_transcript_id,collapsed_transcripts")
    parser.add_argument("output_gtf", help="Output GTF/GFF file path")
    return parser.parse_args()


def load_id_map(map_csv: Path) -> Dict[str, str]:
    id_map: Dict[str, str] = {}
    with map_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"new_transcript_id", "collapsed_transcripts"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError(
                "Map file must contain columns: new_transcript_id, collapsed_transcripts"
            )

        for row in reader:
            new_id = (row.get("new_transcript_id") or "").strip()
            collapsed = (row.get("collapsed_transcripts") or "").strip()
            if not new_id or not collapsed:
                continue

            for old_id in (part.strip() for part in collapsed.split(",")):
                if not old_id:
                    continue
                id_map[old_id] = new_id

    return id_map


def map_id_list(value: str, id_map: Dict[str, str]) -> str:
    items = [item.strip() for item in value.split(",")]
    mapped = [id_map.get(item, item) for item in items]
    return ",".join(mapped)


def remap_gtf_attributes(attr_field: str, id_map: Dict[str, str]) -> str:
    # GTF style: key "value"; key2 "value2";
    if "\"" in attr_field:
        # Keep non-transcript_id attributes byte-for-byte to avoid changing
        # formatting (for example unquoted numeric values).
        segments = re.split(r"(;)", attr_field)
        out_segments: List[str] = []
        for segment in segments:
            if segment == ";":
                out_segments.append(segment)
                continue

            match = re.match(
                r'^(\s*(?:transcript_id|associated_transcript)\s+)(?:"([^"]*)"|(\S+))(\s*)$',
                segment,
            )
            if not match:
                out_segments.append(segment)
                continue

            prefix, quoted_value, unquoted_value, suffix = match.groups()
            original_value = quoted_value if quoted_value is not None else unquoted_value
            mapped_value = id_map.get(original_value, original_value)

            if quoted_value is not None:
                out_segments.append(f'{prefix}"{mapped_value}"{suffix}')
            else:
                out_segments.append(f"{prefix}{mapped_value}{suffix}")

        return "".join(out_segments)

    # GFF3/TLF style: key=value;key2=value2
    parts = [p for p in attr_field.split(";") if p]
    out_parts: List[str] = []
    for part in parts:
        if "=" not in part:
            out_parts.append(part)
            continue
        key, value = part.split("=", 1)
        value = value.strip()
        if key in {"ID", "transcript_id", "associated_transcript", "Parent"}:
            value = map_id_list(value, id_map)
        out_parts.append(f"{key}={value}")
    return ";".join(out_parts)


def remap_file(input_gtf: Path, output_gtf: Path, id_map: Dict[str, str]) -> None:
    with input_gtf.open("r", encoding="utf-8") as in_handle, output_gtf.open(
        "w", encoding="utf-8"
    ) as out_handle:
        for raw_line in in_handle:
            if raw_line.startswith("#") or not raw_line.strip():
                out_handle.write(raw_line)
                continue

            line = raw_line.rstrip("\n")
            fields = line.split("\t")
            if len(fields) < 9:
                out_handle.write(raw_line)
                continue

            fields[8] = remap_gtf_attributes(fields[8], id_map)
            out_handle.write("\t".join(fields) + "\n")


def main() -> int:
    args = parse_args()
    input_gtf = Path(args.input_gtf)
    map_csv = Path(args.map_csv)
    output_gtf = Path(args.output_gtf)

    id_map = load_id_map(map_csv)
    remap_file(input_gtf, output_gtf, id_map)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())