#!/usr/bin/env python3
"""
Test script for the TPBWT cohort splitting CLI tool.
"""

import sys
import os
import tempfile
import shutil

# Add the scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'scripts'))

from tpbwt_cohort_split import get_samples, calculate_optimal_splits


def test_sample_extraction():
    """Test sample extraction from VCF file."""
    print("Testing sample extraction...")
    vcf_path = "tests/data/test.vcf"
    
    try:
        samples = get_samples(vcf_path)
        print(f"✓ Successfully extracted {len(samples)} samples: {samples}")
        return True
    except Exception as e:
        print(f"✗ Failed to extract samples: {e}")
        return False


def test_split_calculation():
    """Test optimal split calculation."""
    print("\nTesting split calculation...")
    
    test_cases = [
        # (num_samples, max_processors, max_memory_per_chunk, expected_result_or_error)
        (100, 16, 10000, 1),  # Should work: processors allow 5, memory needs 1, result is 1
        (50, 32, 5000, 1),    # Should work: processors allow 7, memory needs 1, result is 1  
        (1000, 8, 2500, "error"),  # Should fail: processors allow 3, memory needs 20 - conflict!
        (20, 10, 100, 2),     # Should work: processors allow 4, memory needs 2, result is 2
    ]
    
    all_passed = True
    for num_samples, max_processors, max_memory, expected in test_cases:
        try:
            result = calculate_optimal_splits(num_samples, max_processors, max_memory)
            
            if expected == "error":
                print(f"✗ Expected error for samples={num_samples}, processors={max_processors}, memory={max_memory}, but got result: {result}")
                all_passed = False
            elif result == expected:
                print(f"✓ Samples: {num_samples}, Processors: {max_processors}, Memory: {max_memory} → Splits: {result}")
            else:
                print(f"✗ Samples: {num_samples}, Processors: {max_processors}, Memory: {max_memory} → Splits: {result} (expected {expected})")
                all_passed = False
                
        except ValueError as e:
            if expected == "error":
                print(f"✓ Expected constraint conflict for samples={num_samples}, processors={max_processors}, memory={max_memory}: {str(e)}")
            else:
                print(f"✗ Unexpected error for samples={num_samples}, processors={max_processors}, memory={max_memory}: {e}")
                all_passed = False
        except Exception as e:
            print(f"✗ Failed calculation for samples={num_samples}, processors={max_processors}, memory={max_memory}: {e}")
            all_passed = False
    
    return all_passed


def test_cli_help():
    """Test that the CLI help works."""
    print("\nTesting CLI help...")
    
    import subprocess
    try:
        result = subprocess.run([
            sys.executable, 
            "scripts/tpbwt_cohort_split.py", 
            "--help"
        ], capture_output=True, text=True, cwd=os.path.dirname(__file__))
        
        if result.returncode == 0 and "Split large cohorts" in result.stdout:
            print("✓ CLI help works correctly")
            return True
        else:
            print(f"✗ CLI help failed: {result.stderr}")
            return False
    except Exception as e:
        print(f"✗ CLI help test failed: {e}")
        return False


def main():
    """Run all tests."""
    print("Running TPBWT Cohort Split CLI Tests")
    print("=" * 40)
    
    tests = [
        test_sample_extraction,
        test_split_calculation, 
        test_cli_help
    ]
    
    passed = 0
    for test in tests:
        if test():
            passed += 1
    
    print(f"\nResults: {passed}/{len(tests)} tests passed")
    
    if passed == len(tests):
        print("🎉 All tests passed!")
        return 0
    else:
        print("❌ Some tests failed!")
        return 1


if __name__ == "__main__":
    exit(main())