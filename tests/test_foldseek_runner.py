import unittest
from src.foldseek_runner import check_foldseek_availability, run_foldseek
import os

class TestFoldSeekRunner(unittest.TestCase):
    def test_check_foldseek_availability(self):
        # Assuming foldseek is available on the system
        self.assertEqual(check_foldseek_availability(None), "foldseek")

    def test_run_foldseek(self):
        # This is more of an integration test and may require a mock for subprocess
        # You can use the unittest.mock module to mock subprocess.run
        pass

if __name__ == '__main__':
    unittest.main()
