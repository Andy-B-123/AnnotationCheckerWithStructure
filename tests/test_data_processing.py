import unittest
import pandas as pd
from src.data_processing import protein_recs, filter_genes_transcripts

class TestDataProcessing(unittest.TestCase):
    def test_protein_recs(self):
        # Add tests for protein_recs function

    def test_filter_genes_transcripts(self):
        protein_df = pd.DataFrame({
            'Gene_ID': ['gene1', 'gene1', 'gene2'],
            'Transcript_ID': ['trans1', 'trans2', 'trans1'],
            'sequence': ['MSEQUENCE*', 'MSEQ*', 'MSEQUENCE*'],
            'protein_length': [10, 5, 10]
        })
        filtered_df = filter_genes_transcripts(protein_df, min_length=5, max_length=10, longest_only=True)
        self.assertEqual(len(filtered_df), 2)
        self.assertIn('gene1_trans1', filtered_df['Gene_ID_Transcript_ID'].values)

if __name__ == '__main__':
    unittest.main()
