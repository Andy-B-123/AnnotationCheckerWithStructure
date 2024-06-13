import unittest
from src.io_utils import read_fasta, read_gff
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class TestIOUtils(unittest.TestCase):
    def test_read_fasta(self):
        # Create a temporary fasta file
        test_fasta = "test.fasta"
        records = [SeqRecord(Seq("ATGCGT"), id="test1"), SeqRecord(Seq("ATGCGTAA"), id="test2")]
        with open(test_fasta, "w") as handle:
            SeqIO.write(records, handle, "fasta")

        ref_recs = read_fasta(test_fasta)
        self.assertEqual(len(ref_recs), 2)
        self.assertIn("test1", ref_recs)
        self.assertIn("test2", ref_recs)

        # Clean up
        os.remove(test_fasta)

    # You can add more tests for read_gff, etc.


if __name__ == '__main__':
    unittest.main()
