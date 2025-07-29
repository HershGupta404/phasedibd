#!/usr/bin/env python3
"""
Integration test for the TPBWT cohort splitting CLI tool.
Tests the complete workflow with a small sample dataset.
"""

import sys
import os
import tempfile
import shutil
import subprocess


def create_test_vcf_and_map(temp_dir):
    """Create a small test VCF and genetic map for testing."""
    vcf_path = os.path.join(temp_dir, "test_cohort.vcf")
    map_path = os.path.join(temp_dir, "test_cohort.map")
    
    # Create a VCF with 10 samples
    vcf_content = """##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S001	S002	S003	S004	S005	S006	S007	S008	S009	S010
1	1000	rs1	A	T	999	PASS	.	GT	0|0	0|1	1|0	1|1	0|0	0|1	1|0	1|1	0|0	0|1
1	2000	rs2	G	C	999	PASS	.	GT	1|1	0|0	1|1	1|0	0|0	1|0	0|0	1|1	0|0	1|1
1	3000	rs3	C	T	999	PASS	.	GT	0|1	1|0	0|1	1|1	0|1	1|1	1|0	0|1	1|0	0|1
1	4000	rs4	T	G	999	PASS	.	GT	1|0	1|1	0|0	0|1	1|1	1|1	0|1	1|0	1|1	0|0
1	5000	rs5	A	G	999	PASS	.	GT	0|0	1|0	1|1	1|0	1|1	1|1	0|1	0|0	1|0	1|1
"""
    
    # Create a genetic map
    map_content = """1 rs1 0.1 1000
1 rs2 0.2 2000
1 rs3 0.3 3000
1 rs4 0.4 4000
1 rs5 0.5 5000
"""
    
    with open(vcf_path, 'w') as f:
        f.write(vcf_content)
    
    with open(map_path, 'w') as f:
        f.write(map_content)
    
    return vcf_path, map_path


def test_dry_run():
    """Test the complete workflow without running actual TPBWT analysis."""
    print("Testing complete workflow (dry run)...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test data
        vcf_path, map_path = create_test_vcf_and_map(temp_dir)
        output_dir = os.path.join(temp_dir, "output")
        
        print(f"  - Created test VCF: {vcf_path}")
        print(f"  - Created test map: {map_path}")
        print(f"  - Output directory: {output_dir}")
        
        # Test sample extraction
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'scripts'))
        from tpbwt_cohort_split import get_samples, calculate_optimal_splits
        
        samples = get_samples(vcf_path)
        print(f"  - Extracted {len(samples)} samples: {samples}")
        
        # Test split calculation
        optimal_splits = calculate_optimal_splits(len(samples), 4, 25)  # 10 samples, 4 processors, memory for 25 samples^2
        print(f"  - Calculated {optimal_splits} optimal splits")
        
        # Test CLI argument validation (should not run TPBWT due to missing dependencies)
        cmd = [
            sys.executable,
            "scripts/tpbwt_cohort_split.py",
            "--vcf", vcf_path,
            "--map", map_path,
            "--out", output_dir,
            "--num-splits", "2",
            "--verbose"
        ]
        
        print(f"  - Testing CLI with command: {' '.join(cmd)}")
        
        # We expect this to fail at the TPBWT step, but the preprocessing should work
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.dirname(__file__), timeout=30)
            
            # Check if it got to the point of trying to run TPBWT analysis
            if "Splitting" in result.stderr or "chunks" in result.stderr or "Created filtered VCF" in result.stderr:
                print("  ✓ CLI preprocessing works correctly")
                return True
            elif result.returncode == 0:
                print("  ✓ CLI completed successfully")
                return True
            else:
                print(f"  ✗ CLI failed unexpectedly: {result.stderr}")
                return False
                
        except subprocess.TimeoutExpired:
            print("  ✗ CLI timed out")
            return False
        except Exception as e:
            print(f"  ✗ CLI test failed: {e}")
            return False


def main():
    """Run integration test."""
    print("Running TPBWT Cohort Split CLI Integration Test")
    print("=" * 50)
    
    if test_dry_run():
        print("\n🎉 Integration test passed!")
        return 0
    else:
        print("\n❌ Integration test failed!")
        return 1


if __name__ == "__main__":
    exit(main())