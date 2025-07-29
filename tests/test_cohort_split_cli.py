#!/usr/bin/env python3
"""
Test suite for the TPBWT cohort splitting CLI tool.
"""

import sys
import os
import unittest
import tempfile
import subprocess

# Add the project root to the path
project_root = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, project_root)

from scripts.tpbwt_cohort_split import get_samples, calculate_optimal_splits


class TestCohortSplitCLI(unittest.TestCase):
    """Test cases for the cohort splitting CLI tool."""
    
    def test_sample_extraction(self):
        """Test sample extraction from VCF file."""
        vcf_path = os.path.join(project_root, "tests/data/test.vcf")
        samples = get_samples(vcf_path)
        
        # The test VCF should have 7 samples
        self.assertEqual(len(samples), 7)
        self.assertIn('50', samples)
        self.assertIn('600', samples)
    
    def test_split_calculation_valid(self):
        """Test optimal split calculation for valid constraints."""
        # Test case where memory allows fewer splits than processors
        result = calculate_optimal_splits(100, 16, 10000)
        self.assertEqual(result, 1)
        
        # Test case where both constraints are reasonable
        result = calculate_optimal_splits(20, 10, 100)
        self.assertEqual(result, 2)
    
    def test_split_calculation_conflict(self):
        """Test that conflicting constraints raise ValueError."""
        with self.assertRaises(ValueError) as context:
            calculate_optimal_splits(1000, 8, 2500)
        
        self.assertIn("Cannot satisfy both constraints", str(context.exception))
    
    def test_cli_help(self):
        """Test that the CLI help command works."""
        script_path = os.path.join(project_root, "scripts/tpbwt_cohort_split.py")
        result = subprocess.run([
            sys.executable, script_path, "--help"
        ], capture_output=True, text=True)
        
        self.assertEqual(result.returncode, 0)
        self.assertIn("Split large cohorts", result.stdout)


if __name__ == '__main__':
    unittest.main()