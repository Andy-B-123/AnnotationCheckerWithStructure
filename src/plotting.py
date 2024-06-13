import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

def plot_histogram(protein_df, pdf):
    plt.figure(figsize=(10, 6))
    sns.histplot(protein_df['protein_length'], bins=30, kde=True)
    plt.title('Distribution of Protein Sizes')
    plt.xlabel('Protein Size (amino acids)')
    plt.ylabel('Frequency')
    pdf.savefig()
    plt.close()

def plot_transcripts_per_gene(protein_df, pdf):
    transcripts_per_gene = protein_df.groupby('Gene_ID').size().reset_index(name='transcript_count')
    plt.figure(figsize=(10, 6))
    sns.boxplot(x=transcripts_per_gene['transcript_count'])
    plt.title('Number of Transcripts per Gene')
    plt.xlabel('Number of transcripts per gene')
    pdf.savefig()
    plt.close()

def create_report(protein_df, output_file='report.pdf'):
    with PdfPages(output_file) as pdf:
        plot_histogram(protein_df, pdf)
        plot_transcripts_per_gene(protein_df, pdf)

        valid_count = protein_df['is_valid'].sum()
        invalid_count = len(protein_df) - valid_count
        valid_proportion = valid_count / len(protein_df)
        invalid_proportion = invalid_count / len(protein_df)

        gene_validity = protein_df.groupby('Gene_ID')['is_valid'].any().reset_index()
        gene_invalid_count = len(gene_validity) - gene_validity['is_valid'].sum()
        gene_valid_count = len(gene_validity) - gene_invalid_count
        gene_valid_proportion = gene_valid_count / len(gene_validity)
        gene_invalid_proportion = gene_invalid_count / len(gene_validity)

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
