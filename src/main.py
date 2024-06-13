import argparse
from pathlib import Path
import pandas as pd
from .io_utils import read_fasta
from .data_processing import protein_recs, process_protein_df, filter_genes_transcripts, write_to_fasta
from .plotting import create_report
from .foldseek_runner import run_foldseek

def parse_args():
    parser = argparse.ArgumentParser(description="Generate protein sequences from genome and GFF annotation.")
    parser.add_argument("gff_file", help="Input GFF file with annotations")
    parser.add_argument("ref_file", help="Reference genome file in FASTA format.")
    parser.add_argument("output_directory", help="Provide path to a directory for storing output.")
    parser.add_argument("--foldseek_path", help="Path to the foldseek executable if not available in system path.")
    parser.add_argument("--foldseek_database_path", help="Path to the foldseek database to search")
    parser.add_argument("--threads", type=int, default=1, help="Threads to run FoldSeek")
    return parser.parse_args()

def main():
    args = parse_args()
    ref_recs = read_fasta(args.ref_file)
    protein_df = pd.DataFrame(protein_recs(args.gff_file, ref_recs), columns=['Gene_ID', 'Transcript_ID', 'sequence'])
    protein_df = process_protein_df(protein_df)

    output_directory = Path(args.output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    create_report(protein_df, output_file=output_directory / "InitialAnnotationAssessment.pdf")

    filtered_df = filter_genes_transcripts(protein_df, min_length=5, max_length=1000, longest_only=True)
    output_filtered_fasta_file = output_directory / "filtered_proteins.fasta"
    write_to_fasta(filtered_df, output_filtered_fasta_file)

    if args.foldseek_database_path:
        run_foldseek(output_filtered_fasta_file, args.foldseek_database_path, output_directory, foldseek_path=args.foldseek_path, threads=args.threads)

if __name__ == "__main__":
    main()
