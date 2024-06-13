from __future__ import with_statement
import sys
import os
import operator
from functools import reduce

from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from BCBio import GFF
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path

import subprocess


def main(gff_file, ref_file):
    with open(ref_file) as in_handle:
        ref_recs = SeqIO.to_dict(SeqIO.parse(in_handle, "fasta"))

    protein_df = pd.DataFrame(protein_recs(gff_file, ref_recs), columns=['Gene_ID', 'Transcript_ID', 'sequence'])
    protein_df['protein_length'] = protein_df['sequence'].apply(len)
    protein_df['starts_with_M'] = protein_df['sequence'].apply(lambda seq: seq.startswith('M'))
    protein_df['ends_with_stop'] = protein_df['sequence'].apply(lambda seq: seq.endswith('*'))
    protein_df['no_internal_stops'] = protein_df['sequence'].apply(lambda seq: seq.count('*') == 1)
    protein_df['is_valid'] = protein_df.apply(
        lambda row: row['starts_with_M'] and row['ends_with_stop'] and row['no_internal_stops'], axis=1)

    print(protein_df)
    return protein_df


def protein_recs(gff_file, ref_recs):
    """Generate protein records from GFF gene predictions."""
    with open(gff_file) as in_handle:
        for rec in GFF.parse(in_handle, base_dict=ref_recs):
            for feature in rec.features:
                if feature.type == "gene":
                    gene_id = feature.qualifiers.get("ID", [""])[0]
                    for sub_feature in feature.sub_features:
                        if sub_feature.type == "mRNA":
                            seq_exons = []
                            for cds in sub_feature.sub_features:
                                if cds.type == "CDS" and sub_feature.location.strand != -1:
                                    seq_exons.append(rec.seq[
                                                     int(cds.location.start):
                                                     int(cds.location.end)])
                                elif cds.type == "CDS" and sub_feature.location.strand == -1:
                                    seq_exons.append(rec.seq[
                                                     int(cds.location.start):
                                                     int(cds.location.end)].reverse_complement())
                                else:
                                    continue
                            gene_seq = Seq.Seq(str(reduce(operator.add, seq_exons, "")))
                            protein_seq = gene_seq.translate()
                            transcript_id = sub_feature.qualifiers["ID"][0]
                            yield gene_id, transcript_id, str(protein_seq)


def plot_histogram(protein_df, pdf):
    """Plot a histogram of the distribution of protein sizes."""
    plt.figure(figsize=(10, 6))
    sns.histplot(protein_df['protein_length'], bins=30, kde=True)
    plt.title('Distribution of Protein Sizes')
    plt.xlabel('Protein Size (amino acids)')
    plt.ylabel('Frequency')
    pdf.savefig()
    plt.close()


def plot_transcripts_per_gene(protein_df, pdf):
    """Plot a box-and-whisker plot summarizing the number of transcripts per gene."""
    transcripts_per_gene = protein_df.groupby('Gene_ID').size().reset_index(name='transcript_count')
    plt.figure(figsize=(10, 6))
    sns.boxplot(x=transcripts_per_gene['transcript_count'])
    plt.title('Number of Transcripts per Gene')
    plt.xlabel('Number of transcripts per gene')
    pdf.savefig()
    plt.close()

def create_report(protein_df, output_file='report.pdf'):
    """Create a PDF report with the generated plots and summaries."""
    with PdfPages(output_file) as pdf:
        plot_histogram(protein_df, pdf)
        plot_transcripts_per_gene(protein_df, pdf)

        # Calculate validity for transcripts
        valid_count = protein_df['is_valid'].sum()
        invalid_count = len(protein_df) - valid_count
        valid_proportion = valid_count / len(protein_df)
        invalid_proportion = invalid_count / len(protein_df)

        # Calculate validity for genes
        gene_validity = protein_df.groupby('Gene_ID')['is_valid'].any().reset_index()
        gene_invalid_count = len(gene_validity) - gene_validity['is_valid'].sum()
        gene_valid_count = len(gene_validity) - gene_invalid_count
        gene_valid_proportion = gene_valid_count / len(gene_validity)
        gene_invalid_proportion = gene_invalid_count / len(gene_validity)

        # Write summary to a text page in the PDF
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.axis('off')
        summary_text = (
            f"Count and Proportion of Genes with Only Valid Transcripts:\n"
            f"Count: {gene_valid_count}\n"
            f"Proportion: {gene_valid_proportion:.2%}\n\n"
            f"Count and Proportion of Genes with Any Invalid Transcripts:\n"
            f"Count: {gene_invalid_count}\n"
            f"Proportion: {gene_invalid_proportion:.2%}\n"
            f"Count and Proportion of Valid Transcripts:\n"
            f"Count: {valid_count}\n"
            f"Proportion: {valid_proportion:.2%}\n\n"
            f"Count and Proportion of Invalid Transcripts:\n"
            f"Count: {invalid_count}\n"
            f"Proportion: {invalid_proportion:.2%}\n\n"
        )
        ax.text(0.5, 0.5, summary_text, ha='right', va='center', fontsize=8, wrap=True)
        pdf.savefig()
        plt.close()

    print(f"Initial assessment of annotation saved to {output_file}")

def filter_genes_transcripts(protein_df, min_length, max_length, longest_only):
    """
    Filter the protein_df to keep only transcripts between a minimum and maximum length,
    and filter to keep only the longest if the relevant flag is set to True.

    Args:
    protein_df (pd.DataFrame): DataFrame containing protein data.
    min_length (int): Minimum length of transcripts to keep.
    max_length (int): Maximum length of transcripts to keep.
    longest_only (bool): Flag to keep only the longest transcript per gene.

    Returns:
    pd.DataFrame: Filtered DataFrame.
    """

    # Filter by length
    filtered_df = protein_df[(protein_df['protein_length'] >= min_length) & (protein_df['protein_length'] <= max_length)]

    if longest_only:
        # Keep only the longest transcript per gene
        filtered_df = filtered_df.loc[filtered_df.groupby('Gene_ID')['protein_length'].idxmax()]

    return filtered_df

def write_to_fasta(filtered_df, output_fasta):
    """
    Write the filtered DataFrame to a fasta file with headers as 'Gene_Transcript'.

    Args:
    filtered_df (pd.DataFrame): DataFrame containing filtered protein data.
    output_fasta (str): Output fasta file path.
    """

    records = []
    for _, row in filtered_df.iterrows():
        header = f"{row['Gene_ID']}_{row['Transcript_ID']}"
        sequence = Seq.Seq(row['sequence'])
        record = SeqRecord(sequence, id=header, description="")
        records.append(record)

    with open(output_fasta, "w") as fasta_out:
        SeqIO.write(records, fasta_out, "fasta")

def check_foldseek_availability(foldseek_path):
    """
    Check if foldseek is available on the command line or system path.
    If not, check if the foldseek_path parameter is provided and valid.

    Args:
    foldseek_path (str): Path to the foldseek executable.

    Returns:
    str: Path to the foldseek executable if valid, else raises an error.
    """
    if os.system("command -v foldseek > /dev/null 2>&1") == 0:
        return "foldseek"
    elif foldseek_path and os.path.isfile(foldseek_path) and os.access(foldseek_path, os.X_OK):
        return foldseek_path
    else:
        raise FileNotFoundError("foldseek is not available on the system path or the provided path is invalid.")

def runFoldSeek(fasta_file,database_path, foldseek_path=None):
    foldseek_exec = check_foldseek_availability(foldseek_path)
    # Add code here to run foldseek with the fasta_file and foldseek_exec

    # Define the parameters for the foldseek command
    database_path = database_path
    weights_path = r"/scratch3/bac050/AnnotationCheckerWithStructure/weights"
    output_file = os.path.join(args.output_directory, "filtered_proteins.foldseek.out")
    tmp_dir = os.path.join(args.output_directory, "filtered_proteins.foldseek.tmp")

    # Construct the foldseek command
    command = [
        foldseek_exec, "easy-search", fasta_file,
        database_path, output_file, tmp_dir,
        "--prostt5-model", weights_path, "--format-mode", "2", "--threads", str(threads)
    ]
    print(f"Executing command: {' '.join(command)}")

    # Execute the command
    try:
        subprocess.run(command, check=True)
        print(f"Command executed successfully: {' '.join(command)}")
    except subprocess.CalledProcessError as e:
        print(f"Command failed with error: {e}")


def parse_args():
    parser = argparse.ArgumentParser(description="Generate protein sequences from genome and GFF annotation.")
    parser.add_argument("gff_file", help="Input GFF file with annotations")
    parser.add_argument("ref_file", help="Reference genome file in FASTA format.")
    parser.add_argument("output_directory", help="Provide path to a directory for storing output.")
    parser.add_argument("--foldseek_path", help="Path to the foldseek executable if not available in system path.")
    parser.add_argument("--foldseek_database_path", help="Path to the foldseek database to search")
    parser.add_argument("--threads", help="Threads to run FoldSeek at")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    protein_df = main(args.gff_file, args.ref_file)
#    path_genome = r"U:\Chap4\Sf9_denovo_transcriptome\references\GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.fna"
#    path_annotation = r"U:\Chap4\Sf9_denovo_transcriptome\references\Spod.CSIRO.temp.gff"
#    path_output_directory = r"U:\AnnotationCheckerWithStructure\Development\AnnotationCheckerWithStructure\AnnotationCheckerWithStructure"
    foldseek_path = r"/scratch3/bac050/AnnotationCheckerWithStructure/foldseek/bin/foldseek"
    foldseek_database_path = r"/scratch3/bac050/AnnotationCheckerWithStructure/database_AF-Swiss-prot/AF_Swiss-prot"
    threads = args.threads

    #protein_df = main(path_annotation, path_genome)
    output_directory = Path(args.output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    create_report(protein_df, output_file= os.path.join(args.output_directory,"InitialAnnotationAssessment.pdf"))
#    create_report(protein_df, output_file=os.path.join(path_output_directory, "InitialAnnotationAssessment.pdf"))
    min_length = 5
    max_length = 1000
    longest_only = True
    filtered_df = filter_genes_transcripts(protein_df, min_length, max_length, longest_only)
    output_filtered_fasta_file = os.path.join(args.output_directory,"filtered_proteins.fasta")
#    output_filtered_fasta_file = os.path.join(path_output_directory, "filtered_proteins.fasta")
    write_to_fasta(filtered_df, output_filtered_fasta_file)
    runFoldSeek(output_filtered_fasta_file, foldseek_database_path, foldseek_path=foldseek_path)
    print(os.getcwd())
