#!/usr/bin/env python3
import argparse
import re
import sys
from typing import List, Optional

def parse_args():
    parser = argparse.ArgumentParser(description="Modify AMAS partition file by replacing names with sequence list and making codon partitions based on a pattern")
    parser.add_argument("-i", "--input", required=True, help="Input NEXUS file out of AMAS")
    parser.add_argument("-o", "--output", required=True, help="Output NEXUS file")
    parser.add_argument("-l", "--list", help="List of new partition names in the SAME ORDER as input files given to AMAS")
    parser.add_argument("-c", "--cds", help="Grep-like pattern to identify coding partitions")
    return parser.parse_args()

def sanitize_name(name: str) -> str:
    return re.sub(r"[-./=]", "_", name)

def load_list(list_file: str) -> List[str]:
    with open(list_file) as f:
        return [sanitize_name(line.strip()) for line in f if line.strip()]

def process_nexus(input_file, output_file, name_list=None, cds_pattern=None):
    cds_re = re.compile(cds_pattern) if cds_pattern else None
    new_lines = []
    name_index = 0

    with open(input_file) as f:
        for line in f:
            line_stripped = line.strip()
            if line_stripped.lower().startswith("charset"):
                # Extract name, start, end
                m = re.match(r"charset\s+(\S+)\s*=\s*(\d+)\s*-\s*(\d+);", line_stripped, re.IGNORECASE)
                if not m:
                    new_lines.append(line)
                    continue
                old_name, start, end = m.groups()
                start, end = int(start), int(end)

                # Rename partition
                if name_list:
                    if name_index >= len(name_list):
                        sys.exit("Error: more partitions than lines in --list file")
                    new_name = name_list[name_index]
                else:
                    new_name = sanitize_name(old_name)

                name_index += 1

                # Coding partition
                if cds_re and cds_re.search(new_name):
                    for shift in range(3):
                        new_lines.append(f"        charset {new_name}_c{shift+1} = {start+shift}-{end}\\3;\n")
                else:
                    new_lines.append(f"        charset {new_name} = {start}-{end};\n")
            else:
                new_lines.append(line)

    with open(output_file, "w") as out:
        out.writelines(new_lines)

def main():
    args = parse_args()
    name_list = load_list(args.list) if args.list else None
    process_nexus(args.input, args.output, name_list, args.cds)

if __name__ == "__main__":
    main()

