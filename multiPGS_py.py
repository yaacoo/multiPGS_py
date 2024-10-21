#!/usr/bin/env python3
import gzip
import sys
import os

# Specify the PGS position column name for the desired genome build 
position_column_name = "hm_pos"  # Check the header of the PGS file! This would be either "chr_position" if the PGS is originally in the required genome build, or "hm_pos" if lifted from hg19 to hg38.

def complement(base):
  """Return the complement of a nucleotide base."""
  complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
  return complement.get(base.upper(), base.upper())

def reverse_complement(allele):
  """Return the reverse complement of an allele sequence."""
  return ''.join(complement(base) for base in reversed(allele))

def is_ambiguous_snp(allele1, allele2):
  """Return True if SNP alleles are ambiguous (A/T or G/C), otherwise False."""
  base1 = allele1.upper()
  base2 = allele2.upper()
  if len(base1) != 1 or len(base2) != 1:
      return False  # Not a SNP (length not equal to 1)
  return complement(base1) == base2

def normalize_chromosome(chrom):
  """Normalize chromosome names to allow matching '1' with 'chr1'."""
  chrom = chrom.strip()
  if chrom.startswith('chr'):
      chrom = chrom[3:]
  return chrom

def read_pgs_score_file(pgs_score_file):
  """
  Read the PGS Catalog score file and return a dictionary mapping variants to effect sizes,
  along with the PGS name extracted from the header.
  The key is a tuple: (CHR, POS, REF, ALT)
  """
  pgs_dict = {}
  pgs_name = None
  total_variants = 0  # Total number of variants in PGS file
  with gzip.open(pgs_score_file, 'rt') as f:
      # Read header and extract pgs_name
      for line in f:
          if line.startswith('#'):
              if '=' in line:
                  # Extract pgs_name from header
                  if line.startswith('#pgs_id='):
                      pgs_name = line.strip().split('=', 1)[1]
              continue  # Skip header lines
          if line.strip() == '':
              continue  # Skip empty lines
          # Assume the first non-header line contains column names
          # Read column names to get indices
          header = line.strip().split('\t')
          col_indices = {}
          for idx, col_name in enumerate(header):
              col_indices[col_name] = idx
          # Now process the rest of the file
          break  # Exit the loop to continue processing data
      else:
          print(f"Error: No data found in PGS file {pgs_score_file}.")
          sys.exit(1)
      # Now process the data lines
      for line in f:
          line = line.strip()
          if line == '':
              continue
          tokens = line.strip().split('\t')
          # Extract required columns
          try:
              chr_pgs = tokens[col_indices['chr_name']]
              pos_pgs = tokens[col_indices[position_column_name]]
              ref_pgs = tokens[col_indices['other_allele']]
              alt_pgs = tokens[col_indices['effect_allele']]
              beta = tokens[col_indices['effect_weight']]
          except KeyError as e:
              print(f"Error: Missing required column in PGS file {pgs_score_file}: {e}")
              sys.exit(1)
          except IndexError:
              # Skip this line if any of the required columns are out of index
              continue
          chr_pgs = normalize_chromosome(chr_pgs)
          pos_pgs = int(pos_pgs)
          beta = float(beta)
          key = (chr_pgs, pos_pgs, ref_pgs, alt_pgs)
          pgs_dict[key] = (beta, line)  # Store the full line for unmatched variants
          total_variants += 1
  if pgs_name is None:
      pgs_name = os.path.basename(pgs_score_file)
  return pgs_dict, pgs_name, total_variants

def calculate_prs(vcf_file, pgs_data_list):
  """Calculate the PRS by processing the VCF file line by line, using the DS field, for multiple PGS files."""
  prs_dict = {pgs_name: 0.0 for _, pgs_name, _ in pgs_data_list}
  matched_variants_dict = {pgs_name: set() for _, pgs_name, _ in pgs_data_list}
  total_variants_pgs_dict = {pgs_name: total_variants_pgs for _, pgs_name, total_variants_pgs in pgs_data_list}
  pgs_dicts = {pgs_name: pgs_dict for pgs_dict, pgs_name, _ in pgs_data_list}

  variant_count = 0   

  with gzip.open(vcf_file, 'rt') as f:
      for line in f:
          if line.startswith('#'):
              if line.startswith('#CHROM'):
                  # Extract sample name if needed
                  header_tokens = line.strip().split('\t')
                  if len(header_tokens) < 10:
                      print("Error: VCF does not contain genotype data.")
                      sys.exit(1)
                  sample_name = header_tokens[9]  # Assuming single-sample VCF
              continue  # Skip header lines

          tokens = line.strip().split('\t')
          if len(tokens) < 10:
              continue  # Skip incomplete lines

          chrom_vcf, pos_vcf, _, ref_vcf, alt_vcf, _, _, _, format_field = tokens[:9]
          sample_field = tokens[9]
          pos_vcf = int(pos_vcf)

          # Parse the FORMAT field to find the index of the 'DS' field
          format_keys = format_field.strip().split(':')
          try:
              ds_index = format_keys.index('DS')
          except ValueError:
              # 'DS' not found in FORMAT field
              continue  # Skip this variant
          # Parse the sample field to get the 'DS' value(s)
          sample_values = sample_field.strip().split(':')
          if len(sample_values) <= ds_index:
              continue  # 'DS' value missing
          ds_value_str = sample_values[ds_index]
          if ds_value_str == '.' or ds_value_str == '':
              continue  # Missing dosage value
          # Split DS values for multi-allelic sites
          ds_values_str = ds_value_str.split(',')
          try:
              ds_values = [float(ds) for ds in ds_values_str]
          except ValueError:
              continue  # Invalid DS values, skip variant
          if len(ds_values) != len(alt_vcf.split(',')):
              # The number of DS values should match the number of ALT alleles
              continue  # Skip this variant

          chrom_vcf_norm = normalize_chromosome(chrom_vcf)
          # Possible alleles in VCF
          alleles_vcf = [ref_vcf] + alt_vcf.split(',')
          alt_alleles_vcf = alt_vcf.split(',')
          # Build keys to match with PGS catalogs
          for idx, alt_allele_vcf in enumerate(alt_alleles_vcf):
              ref_allele = ref_vcf.upper()
              alt_allele = alt_allele_vcf.upper()

              # Check if SNP is ambiguous (A/T or G/C) and exclude if so
              if len(ref_allele) == 1 and len(alt_allele) == 1:
                  if is_ambiguous_snp(ref_allele, alt_allele):
                      continue  # Skip ambiguous SNP

              ds_value = ds_values[idx]  # Get the dosage for this alternate allele

              # For each PGS file
              for pgs_name, pgs_dict in pgs_dicts.items():
                  key = (chrom_vcf_norm, pos_vcf, ref_vcf, alt_allele_vcf)
                  key_swap = (chrom_vcf_norm, pos_vcf, alt_allele_vcf, ref_vcf)
                  key_comp = (chrom_vcf_norm, pos_vcf, reverse_complement(ref_vcf), reverse_complement(alt_allele_vcf))
                  key_comp_swap = (chrom_vcf_norm, pos_vcf, reverse_complement(alt_allele_vcf), reverse_complement(ref_vcf))
                  beta = None
                  matched_key = None
                  if key in pgs_dict:
                      beta, _ = pgs_dict[key]
                      matched_key = key
                  elif key_swap in pgs_dict:
                      beta, _ = pgs_dict[key_swap]
                      beta = -beta  # Effect allele is swapped
                      matched_key = key_swap
                  elif key_comp in pgs_dict:
                      beta, _ = pgs_dict[key_comp]
                      matched_key = key_comp
                  elif key_comp_swap in pgs_dict:
                      beta, _ = pgs_dict[key_comp_swap]
                      beta = -beta  # Effect allele is swapped
                      matched_key = key_comp_swap
                  if beta is not None:
                      # Use the DS value as the dosage of effect allele
                      prs_dict[pgs_name] += ds_value * beta
                      matched_variants_dict[pgs_name].add(matched_key)
          # End of processing for the variant
      print(f"Finished processing {variant_count} variants.")

  return prs_dict, matched_variants_dict, total_variants_pgs_dict

def main():
  if len(sys.argv) < 3 or len(sys.argv) > 7:
      print("Usage: python calculate_prs.py <vcf_file.vcf.gz> <pgs_score_file1.txt.gz> [<pgs_score_file2.txt.gz> ... <pgs_score_file5.txt.gz>]")
      sys.exit(1)
  vcf_file = sys.argv[1]
  pgs_score_files = sys.argv[2:]

  if len(pgs_score_files) > 5:
      print("Error: Maximum of 5 PGS score files allowed.")
      sys.exit(1)

  vcf_basename = os.path.basename(vcf_file).split('.')[0]  # Get basename without extension

  print("Reading PGS catalog score files...")
  pgs_data_list = []
  for pgs_score_file in pgs_score_files:
      print(f"Loading {pgs_score_file}...")
      pgs_dict, pgs_name, total_variants_pgs = read_pgs_score_file(pgs_score_file)
      print(f"Loaded {len(pgs_dict)} score entries from {pgs_name}.")
      pgs_data_list.append((pgs_dict, pgs_name, total_variants_pgs))

  print("Calculating PRS from VCF file...")
  prs_dict, matched_variants_dict, total_variants_pgs_dict = calculate_prs(vcf_file, pgs_data_list)

  # Output results for each PGS file
  for pgs_dict, pgs_name, total_variants_pgs in pgs_data_list:
      prs = prs_dict[pgs_name]
      matched_variants = len(matched_variants_dict[pgs_name])
      overlap_percentage = (matched_variants / total_variants_pgs) * 100 if total_variants_pgs > 0 else 0

      print(f"Results for {pgs_name}:")
      print(f"Polygenic Risk Score (PRS): {prs}")
      print(f"Matched Variants: {matched_variants} out of {total_variants_pgs} ({overlap_percentage:.2f}%)")

      # Output the result to a single-line text file
      output_filename = f"{vcf_basename}_{pgs_name}_PRS_result.txt"
      with open(output_filename, 'w') as out_file:
          out_file.write(f"{vcf_basename}\t{pgs_name}\t{prs}\t{matched_variants}\t{total_variants_pgs}\n")
      print(f"Results written to {output_filename}")

if __name__ == "__main__":
  main()
