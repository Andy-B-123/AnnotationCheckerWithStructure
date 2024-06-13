from functools import reduce
import operator
from Bio import Seq
import pandas as pd
from BCBio import GFF

def protein_recs(gff_file, ref_recs):
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
                                    seq_exons.append(rec.seq[int(cds.location.start): int(cds.location.end)])
                                elif cds.type == "CDS" and sub_feature.location.strand == -1:
                                    seq_exons.append(rec.seq[int(cds.location.start): int(cds.location.end)].reverse_complement())
                            gene_seq = Seq.Seq(str(reduce(operator.add, seq_exons, "")))
                            protein_seq = gene_seq.translate()
                            transcript_id = sub_feature.qualifiers["ID"][0]
                            yield gene_id, transcript_id, str(protein_seq)

def process_protein_df(protein_df):
    protein_df['protein_length'] = protein_df['sequence'].apply(len)
    protein_df['starts_with_M'] = protein_df['sequence'].apply(lambda seq: seq.startswith('M'))
    protein_df['ends_with_stop'] = protein_df['sequence'].apply(lambda seq: seq.endswith('*'))
    protein_df['no_internal_stops'] = protein_df['sequence'].apply(lambda seq: seq.count('*') == 1)
    protein_df['is_valid'] = protein_df.apply(lambda row: row['starts_with_M'] and row['ends_with_stop'] and row['no_internal_stops'], axis=1)
    return protein_df

def filter_genes_transcripts(protein_df, min_length, max_length, longest_only):
    filtered_df = protein_df[(protein_df['protein_length'] >= min_length) & (protein_df['protein_length'] <= max_length)]
    if longest_only:
        filtered_df = filtered_df.loc[filtered_df.groupby('Gene_ID')['protein_length'].idxmax()]
    return filtered_df

def write_to_fasta(filtered_df, output_fasta):
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    from Bio.Seq import Seq
    records = [SeqRecord(Seq(row['sequence']), id=f"{row['Gene_ID']}_{row['Transcript_ID']}", description="") for _, row in filtered_df.iterrows()]
    with open(output_fasta, "w") as fasta_out:
        SeqIO.write(records, fasta_out, "fasta")
