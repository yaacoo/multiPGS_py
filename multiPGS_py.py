#!/usr/bin/env python3
import pysam
import gzip
import sys
import os
from typing import Dict, Set, Tuple
import logging
from collections import defaultdict

position_column_name = "hm_pos" # Please check the PGS file for the correct column name that match your VCF's genome build.

# Set up logging
logging.basicConfig(
  level=logging.INFO,
  format='%(asctime)s - %(levelname)s - %(message)s'
)

def complement(base: str) -> str:
  """Return the complement of a nucleotide base."""
  COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
  return COMPLEMENT.get(base.upper(), base.upper())

def reverse_complement(allele: str) -> str:
  """Return the reverse complement of an allele sequence."""
  return ''.join(complement(base) for base in reversed(allele))

def is_ambiguous_snp(allele1: str, allele2: str) -> bool:
  """Return True if SNP alleles are ambiguous (A/T or G/C)."""
  if len(allele1) != 1 or len(allele2) != 1:
      return False
  return complement(allele1.upper()) == allele2.upper()

def normalize_chromosome(chrom: str) -> str:
  """Normalize chromosome names."""
  return chrom[3:] if chrom.startswith('chr') else chrom

def read_pgs_score_file(pgs_score_file: str) -> Tuple[Dict, str, int]:
  """Read PGS file and return dictionary of variants with their effect sizes."""
  pgs_dict = {}
  pgs_name = None
  total_variants = 0
  
  with gzip.open(pgs_score_file, 'rt') as f:
      # Process header
      for line in f:
          if line.startswith('#'):
              if line.startswith('#pgs_id='):
                  pgs_name = line.strip().split('=', 1)[1]
              continue
          if not line.strip():
              continue
          
          header = line.strip().split('\t')
          col_indices = {col: idx for idx, col in enumerate(header)}
          break
      
      # Process variants
      for line in f:
          if not line.strip():
              continue
          
          tokens = line.strip().split('\t')
          try:
              chr_pgs = normalize_chromosome(tokens[col_indices['chr_name']])
              pos_pgs = int(tokens[col_indices[position_column_name]])
              ref_pgs = tokens[col_indices['other_allele']]
              alt_pgs = tokens[col_indices['effect_allele']]
              beta = float(tokens[col_indices['effect_weight']])
              
              # Store variant information
              key = (chr_pgs, pos_pgs, ref_pgs, alt_pgs)
              pgs_dict[key] = beta
              total_variants += 1
              
          except (KeyError, IndexError, ValueError) as e:
              continue
  
  return pgs_dict, pgs_name or os.path.basename(pgs_score_file), total_variants

def create_position_lookup(pgs_dict: Dict) -> Dict[str, Set[int]]:
  """Create a chromosome-based lookup of positions needed."""
  positions_by_chrom = defaultdict(set)
  for chrom, pos, _, _ in pgs_dict.keys():
      positions_by_chrom[chrom].add(pos)
  return positions_by_chrom

def process_variant(record, pgs_dict: Dict, pgs_name: str, 
                 prs_score: float, matched_variants: Set) -> float:
  """Process a single variant and update PRS score."""
  chrom = normalize_chromosome(record.chrom)
  pos = record.pos
  ref = record.ref
  
  # Get dosage value
  if 'DS' not in record.format:
      return 0.0
  
  ds_values = record.samples[0]['DS']
  if not isinstance(ds_values, (list, tuple)):
      ds_values = [ds_values]
  
  score_contribution = 0.0
  
  # Process each alternate allele
  for alt_idx, alt in enumerate(record.alts):
      if alt_idx >= len(ds_values):
          continue
          
      ds_value = ds_values[alt_idx]
      if ds_value is None:
          continue
          
      ref_allele = ref.upper()
      alt_allele = alt.upper()
      
      # Skip ambiguous SNPs
      if len(ref_allele) == 1 and len(alt_allele) == 1:
          if is_ambiguous_snp(ref_allele, alt_allele):
              continue
      
      # Check all possible variant configurations
      keys = [
          (chrom, pos, ref, alt),
          (chrom, pos, alt, ref),
          (chrom, pos, reverse_complement(ref), reverse_complement(alt)),
          (chrom, pos, reverse_complement(alt), reverse_complement(ref))
      ]
      
      for key in keys:
          if key in pgs_dict:
              beta = pgs_dict[key]
              if key[2] != ref:  # If alleles are swapped
                  beta = -beta
              score_contribution += ds_value * beta
              matched_variants.add(key)
              break
  
  return score_contribution

def calculate_prs_with_tabix(vcf_file: str, pgs_data_list: list) -> Tuple[Dict, Dict, Dict]:
  """Calculate PRS using tabix-based position lookup."""
  prs_dict = {pgs_name: 0.0 for _, pgs_name, _ in pgs_data_list}
  matched_variants_dict = {pgs_name: set() for _, pgs_name, _ in pgs_data_list}
  total_variants_dict = {pgs_name: total for _, pgs_name, total in pgs_data_list}
  
  # Create position lookups for each PGS
  position_lookups = {}
  pgs_dicts = {}
  for pgs_dict, pgs_name, _ in pgs_data_list:
      position_lookups[pgs_name] = create_position_lookup(pgs_dict)
      pgs_dicts[pgs_name] = pgs_dict
  
  try:
      # Open VCF file with tabix
      with pysam.VariantFile(vcf_file) as vcf:
          # Process each chromosome
          for chrom in vcf.header.contigs:
              chrom_norm = normalize_chromosome(chrom)
              logging.info(f"Processing chromosome {chrom}")
              
              # Collect positions needed for this chromosome
              for pgs_name, positions in position_lookups.items():
                  if chrom_norm not in positions:
                      continue
                  
                  # Fetch and process variants
                  for pos in sorted(positions[chrom_norm]):
                      try:
                          for record in vcf.fetch(chrom, pos-1, pos):
                              score = process_variant(
                                  record,
                                  pgs_dicts[pgs_name],
                                  pgs_name,
                                  prs_dict[pgs_name],
                                  matched_variants_dict[pgs_name]
                              )
                              prs_dict[pgs_name] += score
                      except ValueError:
                          continue
                      
  except Exception as e:
      logging.error(f"Error processing VCF file: {e}")
      raise
  
  return prs_dict, matched_variants_dict, total_variants_dict

def main():
  if len(sys.argv) < 3 or len(sys.argv) > 7:
      logging.error("Usage: python calculate_prs.py <vcf_file.vcf.gz> <pgs_score_file1.txt.gz> [<pgs_score_file2.txt.gz> ... <pgs_score_file5.txt.gz>]")
      sys.exit(1)
  
  vcf_file = sys.argv[1]
  pgs_score_files = sys.argv[2:]
  vcf_basename = os.path.basename(vcf_file).split('.')[0]
  
  # Check if VCF file is indexed
  if not os.path.exists(vcf_file + '.tbi'):
      logging.error(f"VCF file {vcf_file} must be bgzipped and indexed with tabix")
      logging.error("Run: tabix -p vcf your_file.vcf.gz")
      sys.exit(1)
  
  # Load PGS files
  logging.info("Loading PGS files...")
  pgs_data_list = []
  for pgs_file in pgs_score_files:
      logging.info(f"Processing {pgs_file}")
      pgs_dict, pgs_name, total_variants = read_pgs_score_file(pgs_file)
      pgs_data_list.append((pgs_dict, pgs_name, total_variants))
      logging.info(f"Loaded {total_variants} variants from {pgs_name}")
  
  # Calculate PRS
  logging.info("Calculating PRS...")
  prs_dict, matched_variants_dict, total_variants_dict = calculate_prs_with_tabix(
      vcf_file, pgs_data_list
  )
  
  # Output results
  for _, pgs_name, _ in pgs_data_list:
      prs = prs_dict[pgs_name]
      matched_variants = len(matched_variants_dict[pgs_name])
      total_variants = total_variants_dict[pgs_name]
      overlap_percentage = (matched_variants / total_variants * 100) if total_variants > 0 else 0
      
      logging.info(f"\nResults for {pgs_name}:")
      logging.info(f"PRS: {prs}")
      logging.info(f"Matched Variants: {matched_variants}/{total_variants} ({overlap_percentage:.2f}%)")
      
      # Write results to file
      output_filename = f"{vcf_basename}_{pgs_name}_PRS_result.txt"
      with open(output_filename, 'w') as out_file:
          out_file.write(f"{vcf_basename}\t{pgs_name}\t{prs}\t{matched_variants}\t{total_variants}\n")
      logging.info(f"Results written to {output_filename}")

if __name__ == "__main__":
  main()
