"""
final_assembly_assessment.py
Marcus Viscardi, 11/6/2024

Quick script to take the contig_counts.txt from samtools idxstats
and use those to pick the contig that was actually sequenced significantly,
based on reads per kb of contig.
"""

import pandas as pd
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("contig_count_file",
                    help="output of samtools idxstats")
parser.add_argument("medaka_assembly_file",
                    help="the medaka polished fasta assembly file")
parser.add_argument("mapped_to_medaka",
                    help="output path to the sorted  bam file of reads mapped to the medaka assembly")
parser.add_argument("final_assembly",
                    help="output path to the chosen assembly contig")
parser.add_argument("assembly_mapped_reads",
                    help="output path to a summary file")
parser.add_argument("--min_reads_per_kb", "-m", type=float, default=10,)

args = parser.parse_args()

if __name__ == '__main__':
    contig_counts = pd.read_table(args.contig_count_file, names=["contig", "seq_len", "map_count", "unmap_count"])

    contig_counts['reads_per_kb'] = contig_counts['map_count'] / contig_counts['seq_len'] * 1000

    contig_counts.sort_values("reads_per_kb", ascending=False, inplace=True)
    contig_counts.drop(columns=["unmap_count"], inplace=True)
    contig_counts['will_be_dropped_based_on_reads_per_kb'] = contig_counts['reads_per_kb'] < args.min_reads_per_kb

    contig_counts.to_csv(args.assembly_mapped_reads, sep="\t", index=False)

    medaka_assembly_df = pd.DataFrame()
    # This assumes the incoming medaka assembly file is in 2 line fasta format:
    with open(args.medaka_assembly_file, 'r') as input_fasta:
        lines = input_fasta.readlines()
        for i in range(0, len(lines), 2):
            contig_name = lines[i].strip()[1:]
            contig_seq = lines[i+1].strip()
            medaka_assembly_df = medaka_assembly_df.append({"contig": contig_name,
                                                            "seq": contig_seq},
                                                           ignore_index=True)
    print(medaka_assembly_df)
    merge_df = pd.merge(contig_counts, medaka_assembly_df, on="contig")
    merge_df = merge_df[merge_df['reads_per_kb'] >= args.min_reads_per_kb]
    merge_df.sort_values("reads_per_kb", ascending=False, inplace=True)

    with open(args.final_assembly, 'w') as output_fasta:
        for i, row in merge_df.iterrows():
            new_tagline = (f">{row['contig']} "
                           f"length={row['seq_len']} "
                           f"read_count={row['map_count']} "
                           f"reads_per_kb={row['reads_per_kb']:0.2f}\n")
            output_fasta.write(new_tagline + row['seq'] + "\n")
