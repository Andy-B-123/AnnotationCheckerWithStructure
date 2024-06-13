import unittest
import pandas as pd
from src.plotting import plot_histogram, plot_transcripts_per_gene
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import tempfile
import os

class TestPlotting(unittest.TestCase):
    def setUp(self):
        self.protein_df = pd.DataFrame({
            'Gene_ID': ['gene1', 'gene1', 'gene2'],
            'Transcript_ID': ['trans1', 'trans2', 'trans1'],
            'sequence': ['MSEQUENCE*', 'MSEQ*', 'MSEQUENCE*'],
            'protein_length': [10, 5, 10]
        })
        self.pdf_file = tempfile.NamedTemporaryFile(delete=False)

    def tearDown(self):
        os.remove(self.pdf_file.name)

    def test_plot_histogram(self):
        with PdfPages(self.pdf_file.name) as pdf:
            plot_histogram(self.protein_df, pdf)
        self.assertTrue(os.path.isfile(self.pdf_file.name))

    def test_plot_transcripts_per_gene(self):
        with PdfPages(self.pdf_file.name) as pdf:
            plot_transcripts_per_gene(self.protein_df, pdf)
        self.assertTrue(os.path.isfile(self.pdf_file.name))

if __name__ == '__main__':
    unittest.main()
