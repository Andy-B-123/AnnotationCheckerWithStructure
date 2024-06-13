from Bio import SeqIO

def read_fasta(ref_file):
    with open(ref_file) as in_handle:
        return SeqIO.to_dict(SeqIO.parse(in_handle, "fasta"))

def read_gff(gff_file):
    with open(gff_file) as in_handle:
        return GFF.parse(in_handle)
