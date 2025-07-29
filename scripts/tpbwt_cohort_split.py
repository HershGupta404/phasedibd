#!/usr/bin/env python3
"""
TPBWT Cohort Splitting Tool

This CLI tool splits large cohorts into smaller groups for efficient IBD computation
based on memory and processor constraints. It combines in-sample and out-of-sample
TPBWT analysis with automatic determination of optimal split counts.

The tool implements the following constraints:
- Given P processors, max splits N: (N+1)(N)/2 = P (out-of-sample comparisons)
- Memory constraint: (S/N)^2 = memory, where S is samples and N is splits
"""

import argparse
import os
import subprocess
import sys
from itertools import combinations
from math import ceil, sqrt
from concurrent.futures import ProcessPoolExecutor
import logging


def configure_logging(verbose=False):
    """Configure logging based on verbosity level."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.StreamHandler()]
    )


def calculate_optimal_splits(num_samples, max_processors, max_memory_per_chunk):
    """
    Calculate the optimal number of splits based on processor and memory constraints.
    
    Args:
        num_samples: Total number of samples in the cohort
        max_processors: Maximum number of available processors
        max_memory_per_chunk: Maximum memory (in samples^2) per chunk
    
    Returns:
        Optimal number of splits
    
    Raises:
        ValueError: If constraints cannot be satisfied
    """
    # Constraint 1: Processor limit - (N+1)(N)/2 <= P
    # Solve N^2 + N - 2P <= 0 using quadratic formula
    # N <= (-1 + sqrt(1 + 8P)) / 2
    max_splits_by_processors = int((-1 + sqrt(1 + 8 * max_processors)) / 2)
    
    # Constraint 2: Memory limit - (S/N)^2 <= memory
    # S^2/N^2 <= memory, so N >= S/sqrt(memory)
    min_splits_by_memory = ceil(num_samples / sqrt(max_memory_per_chunk))
    
    # Check if constraints can be satisfied
    if min_splits_by_memory > max_splits_by_processors:
        raise ValueError(
            f"Cannot satisfy both constraints: "
            f"Memory requires at least {min_splits_by_memory} splits, "
            f"but processors allow at most {max_splits_by_processors} splits. "
            f"Either increase --max-processors or decrease --max-memory."
        )
    
    # Choose the more restrictive constraint (higher minimum, lower maximum)
    optimal_splits = max(min_splits_by_memory, 1)  # At least 1 split
    optimal_splits = min(optimal_splits, max_splits_by_processors)  # But not more than processor limit
    optimal_splits = min(optimal_splits, num_samples)  # And not more than samples
    
    # Ensure we have at least 1 split
    if optimal_splits < 1:
        optimal_splits = 1
    
    logging.info(f"Calculated optimal splits: {optimal_splits}")
    logging.info(f"  - Max splits by processors ({max_processors}): {max_splits_by_processors}")
    logging.info(f"  - Min splits by memory ({max_memory_per_chunk}): {min_splits_by_memory}")
    logging.info(f"  - Out-of-sample comparisons: {(optimal_splits + 1) * optimal_splits // 2}")
    logging.info(f"  - Samples per chunk: ~{ceil(num_samples / optimal_splits)}")
    
    return optimal_splits


def get_samples(vcf_path):
    """Extract sample list from VCF file using plink2 or direct VCF parsing."""
    logging.info(f"Extracting sample list from {vcf_path}")
    
    # Try plink2 first if available
    try:
        cmd = [
            "plink2",
            "--vcf", vcf_path,
            "--write-samples",
            "--out", "temp_samples",
            "--threads", "32"
        ]
        subprocess.run(cmd, check=True)
        samples = []
        with open("temp_samples.id") as f:
            for line in f:
                if line.startswith("FID") or not line.strip():
                    continue
                samples.append(line.strip())
        os.remove("temp_samples.id")
        logging.info(f"Found {len(samples)} samples using plink2")
        return samples
    
    except (subprocess.CalledProcessError, FileNotFoundError):
        logging.info("plink2 not available, parsing VCF header directly")
        
        # Fallback: parse VCF header directly
        samples = []
        with open(vcf_path, 'r') as f:
            for line in f:
                if line.startswith('#CHROM'):
                    # Parse sample names from header line
                    parts = line.strip().split('\t')
                    if len(parts) > 9:  # Standard VCF has 9 fixed columns
                        samples = parts[9:]  # Sample names start from column 10
                    break
                elif not line.startswith('#'):
                    # No header found
                    break
        
        if not samples:
            raise ValueError(f"Could not extract samples from VCF file: {vcf_path}")
        
        logging.info(f"Found {len(samples)} samples using direct VCF parsing")
        return samples


def write_sample_file(samples, path):
    """Write list of samples to a file."""
    with open(path, 'w') as f:
        for s in samples:
            f.write(s + "\n")


def split_vcf_plink(vcf_path, samples, out_prefix):
    """Split VCF file to include only specified samples using plink2 or custom implementation."""
    
    # Try plink2 first if available
    try:
        keep_file = f"{out_prefix}.keep.txt"
        write_sample_file(samples, keep_file)
        logging.info(f"Splitting VCF to {out_prefix}.vcf with {len(samples)} samples using plink2")
        cmd = [
            "plink2",
            "--vcf", vcf_path,
            "--keep", keep_file,
            "--recode", "vcf",
            "--out", out_prefix,
            "--threads", "32"
        ]
        subprocess.run(cmd, check=True)
        os.remove(keep_file)
        return f"{out_prefix}.vcf"
    
    except (subprocess.CalledProcessError, FileNotFoundError):
        logging.info("plink2 not available, using custom VCF splitting")
        
        # Fallback: custom VCF splitting implementation
        return split_vcf_custom(vcf_path, samples, out_prefix)


def split_vcf_custom(vcf_path, samples, out_prefix):
    """Custom VCF splitting implementation when plink2 is not available."""
    out_vcf_path = f"{out_prefix}.vcf"
    samples_set = set(samples)
    
    with open(vcf_path, 'r') as infile, open(out_vcf_path, 'w') as outfile:
        sample_indices = []
        
        for line in infile:
            if line.startswith('##'):
                # Copy header lines
                outfile.write(line)
            elif line.startswith('#CHROM'):
                # Parse and filter sample columns
                parts = line.strip().split('\t')
                fixed_cols = parts[:9]  # First 9 columns are fixed
                all_samples = parts[9:] if len(parts) > 9 else []
                
                # Find indices of samples to keep
                sample_indices = [i for i, sample in enumerate(all_samples) 
                                if sample in samples_set]
                
                # Write filtered header
                kept_samples = [all_samples[i] for i in sample_indices]
                outfile.write('\t'.join(fixed_cols + kept_samples) + '\n')
                
                logging.info(f"Keeping {len(kept_samples)} out of {len(all_samples)} samples")
                
            elif not line.startswith('#'):
                # Filter data lines
                parts = line.strip().split('\t')
                fixed_cols = parts[:9]
                sample_cols = parts[9:] if len(parts) > 9 else []
                
                # Keep only selected sample columns
                kept_sample_cols = [sample_cols[i] for i in sample_indices 
                                  if i < len(sample_cols)]
                
                outfile.write('\t'.join(fixed_cols + kept_sample_cols) + '\n')
    
    logging.info(f"Created filtered VCF: {out_vcf_path}")
    return out_vcf_path


def run_insample_driver(vcf_path, map_path, prefix):
    """Run in-sample TPBWT analysis."""
    logging.info(f"Running in-sample TPBWT: {prefix}")
    
    # Import phasedibd for the in-sample analysis
    try:
        import phasedibd as ibd
        
        haplotypes = ibd.VcfHaplotypeAlignment(vcf_path, map_path)
        tpbwt = ibd.TPBWTAnalysis()
        tpbwt.compute_ibd(
            haplotypes,
            segments_out_path=f"{prefix}_segs",
            compressed_out_path=f"{prefix}_tpbwt"
        )
        logging.info(f"Completed in-sample TPBWT: {prefix}")
    except Exception as e:
        logging.error(f"Error in in-sample analysis for {prefix}: {e}")
        raise


def run_outsample_driver(tpbwt1, tpbwt2, prefix):
    """Run out-of-sample TPBWT analysis."""
    logging.info(f"Running out-sample TPBWT: {prefix}")
    
    def normalize_dir(path):
        return os.path.join(path, '') if os.path.isdir(path) else path
    
    try:
        import phasedibd as ibd
        
        h1 = ibd.CompressedHaplotypeAlignment(normalize_dir(tpbwt1))
        h2 = ibd.CompressedHaplotypeAlignment(normalize_dir(tpbwt2))
        tpbwt = ibd.TPBWTAnalysis()
        tpbwt.compute_ibd(
            h1, h2,
            segments_out_path=f"{prefix}_segs"
        )
        logging.info(f"Completed out-sample TPBWT: {prefix}")
    except Exception as e:
        logging.error(f"Error in out-sample analysis for {prefix}: {e}")
        raise


def create_argument_parser():
    """Create and configure argument parser."""
    parser = argparse.ArgumentParser(
        description="Split large cohorts for efficient TPBWT IBD computation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-determine splits based on constraints
  python tpbwt_cohort_split.py --vcf data.vcf --map genetic.map --out results/ --max-processors 16 --max-memory 10000

  # Manually specify number of splits
  python tpbwt_cohort_split.py --vcf data.vcf --map genetic.map --out results/ --num-splits 4

  # Use all available cores
  python tpbwt_cohort_split.py --vcf data.vcf --map genetic.map --out results/ --max-processors auto --max-memory 10000
        """
    )
    
    # Required arguments
    parser.add_argument("--vcf", required=True, help="Input phased VCF file")
    parser.add_argument("--map", required=True, help="Genetic map file")
    parser.add_argument("--out", required=True, help="Output directory prefix")
    
    # Split configuration
    split_group = parser.add_mutually_exclusive_group(required=True)
    split_group.add_argument("--num-splits", type=int, 
                           help="Manual number of splits (overrides automatic calculation)")
    split_group.add_argument("--max-processors", type=str,
                           help="Maximum processors available ('auto' or integer)")
    
    parser.add_argument("--max-memory", type=int, default=10000,
                       help="Maximum memory per chunk in samples^2 (default: 10000)")
    
    # Processing options
    parser.add_argument("--max-workers", type=int, default=None,
                       help="Maximum parallel workers (default: number of splits)")
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="Enable verbose logging")
    
    return parser


def main():
    """Main execution function."""
    parser = create_argument_parser()
    args = parser.parse_args()
    
    configure_logging(args.verbose)
    
    # Validate and create output directory
    os.makedirs(args.out, exist_ok=True)
    logging.info(f"Output directory: {args.out}")
    
    try:
        # Get samples from VCF
        samples = get_samples(args.vcf)
        num_samples = len(samples)
        
        # Determine number of splits
        if args.num_splits:
            num_splits = args.num_splits
            logging.info(f"Using manual split count: {num_splits}")
        else:
            # Parse max_processors
            if args.max_processors.lower() == 'auto':
                max_processors = os.cpu_count()
            else:
                max_processors = int(args.max_processors)
            
            num_splits = calculate_optimal_splits(num_samples, max_processors, args.max_memory)
        
        # Validate splits
        if num_splits < 1:
            raise ValueError("Number of splits must be at least 1")
        if num_splits > num_samples:
            raise ValueError(f"Number of splits ({num_splits}) cannot exceed number of samples ({num_samples})")
        
        # Split samples into chunks
        chunk_size = ceil(num_samples / num_splits)
        chunks = [samples[i * chunk_size : min((i + 1) * chunk_size, num_samples)] 
                 for i in range(num_splits)]
        
        logging.info(f"Splitting {num_samples} samples into {len(chunks)} chunks")
        for i, chunk in enumerate(chunks, 1):
            logging.info(f"  Chunk {i}: {len(chunk)} samples")
        
        # Create VCF files for each chunk and prepare for processing
        chunk_vcfs = []
        for i, chunk in enumerate(chunks, start=1):
            prefix = os.path.join(args.out, f"chunk{i}")
            vcf_path = split_vcf_plink(args.vcf, chunk, prefix)
            chunk_vcfs.append((vcf_path, args.map, prefix))
        
        # Determine number of workers
        max_workers = args.max_workers if args.max_workers else len(chunks)
        
        # Run in-sample TPBWT jobs
        logging.info("Starting in-sample TPBWT jobs...")
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(run_insample_driver, vcf, mapf, prefix)
                for vcf, mapf, prefix in chunk_vcfs
            ]
            for future in futures:
                future.result()
        logging.info("Completed all in-sample TPBWT jobs")
        
        # Run out-of-sample TPBWT jobs
        compressed_paths = [f"{prefix}_tpbwt" for _, _, prefix in chunk_vcfs]
        logging.info("Starting out-sample TPBWT jobs...")
        
        out_sample_pairs = list(combinations(enumerate(compressed_paths, start=1), 2))
        logging.info(f"Running {len(out_sample_pairs)} out-of-sample comparisons")
        
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            for (i, p1), (j, p2) in out_sample_pairs:
                out_prefix = os.path.join(args.out, f"chunk{i}_chunk{j}")
                futures.append(executor.submit(run_outsample_driver, p1, p2, out_prefix))
            
            for future in futures:
                future.result()
        
        logging.info("Completed all out-sample TPBWT jobs")
        logging.info(f"Results written to: {args.out}")
        
    except subprocess.CalledProcessError as e:
        logging.error(f"Subprocess failed: {e}")
        sys.exit(1)
    except Exception as e:
        logging.exception("Unexpected error occurred")
        sys.exit(1)


if __name__ == "__main__":
    main()