#!/usr/bin/env python3
import os
import sys
import logging
import multiprocessing as mp
from typing import List
from pathlib import Path
import subprocess
from datetime import datetime

n_cpus=20 # Change based on the number of CPUs available, if submitted on a cluster or VM, recommended to allocate 1GB of RAM per CPU 

# Set up logging
logging.basicConfig(
  level=logging.INFO,
  format='%(asctime)s - %(levelname)s - %(message)s',
  handlers=[
      logging.FileHandler('parallel_prs.log'),
      logging.StreamHandler()
  ]
)

def process_vcf(args: tuple) -> str:
  """
  Process a single VCF file with the PRS calculator script
  """
  vcf_path, pgs_files = args
  try:
      cmd = ["python", "multiPGS_py.py", vcf_path] + pgs_files
      result = subprocess.run(
          cmd,
          stdout=subprocess.PIPE,
          stderr=subprocess.PIPE,
          universal_newlines=True
      )
      if result.returncode == 0:
          logging.info(f"Successfully processed {vcf_path}")
          return f"Success: {vcf_path}"
      else:
          logging.error(f"Error processing {vcf_path}: {result.stderr}")
          return f"Failed: {vcf_path} - {result.stderr}"
  except subprocess.CalledProcessError as e:
      logging.error(f"Error processing {vcf_path}: {str(e)}")
      return f"Failed: {vcf_path} - {str(e)}"
  except Exception as e:
      logging.error(f"Unexpected error processing {vcf_path}: {str(e)}")
      return f"Failed: {vcf_path} - {str(e)}"

def process_vcf_batch(vcf_list: List[str], pgs_files: List[str], num_cores: int = n_cpus):
  """
  Process a batch of VCF files in parallel
  """
  with mp.Pool(processes=num_cores) as pool:
      # Create list of arguments for each process
      process_args = [(vcf, pgs_files) for vcf in vcf_list]
      # Map the processing function to the arguments
      results = pool.map(process_vcf, process_args)
  return results

def main():
  if len(sys.argv) < 3:
      print("Usage: python parallel_pgs.py <vcf_list.txt> <pgs_file1.txt.gz> [<pgs_file2.txt.gz> ...]")
      sys.exit(1)

  # Read VCF list file
  vcf_list_file = sys.argv[1]
  pgs_files = sys.argv[2:]

  # Verify PGS files exist
  for pgs_file in pgs_files:
      if not os.path.exists(pgs_file):
          logging.error(f"PGS file not found: {pgs_file}")
          sys.exit(1)

  # Read and validate VCF paths
  try:
      with open(vcf_list_file) as f:
          vcf_paths = [line.strip() for line in f if line.strip()]
  except Exception as e:
      logging.error(f"Error reading VCF list file: {e}")
      sys.exit(1)

  # Validate VCF paths
  valid_vcf_paths = []
  for vcf_path in vcf_paths:
      if os.path.exists(vcf_path):
          if os.path.exists(vcf_path + '.tbi'):
              valid_vcf_paths.append(vcf_path)
          else:
              logging.warning(f"Skipping {vcf_path}: No tabix index found")
      else:
          logging.warning(f"Skipping {vcf_path}: File not found")

  if not valid_vcf_paths:
      logging.error("No valid VCF files found")
      sys.exit(1)

  # Process VCFs in batches
  batch_size = n_cpus
  num_cores = min(n_cpus, len(valid_vcf_paths))  # Use fewer cores if fewer VCFs
  
  logging.info(f"Starting parallel processing of {len(valid_vcf_paths)} VCF files")
  logging.info(f"Using {num_cores} cores")
  
  # Process in batches
  for i in range(0, len(valid_vcf_paths), batch_size):
      batch = valid_vcf_paths[i:i + batch_size]
      logging.info(f"Processing batch {i//batch_size + 1} ({len(batch)} files)")
      
      results = process_vcf_batch(batch, pgs_files, num_cores)
      
      # Log results for this batch
      for result in results:
          logging.info(result)

  logging.info("All VCF files have been processed")

if __name__ == "__main__":
  start_time = datetime.now()
  main()
  end_time = datetime.now()
  logging.info(f"Total execution time: {end_time - start_time}")
