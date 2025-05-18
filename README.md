# multiPGS_py: Fast, simple and low-memory PGS scoring

multiPGS_py is a fast, simple, and low-memory Python method to calculate polygenic scores (PGS/PRS) from any PGS Catalog weights (pgscatalog.org) and single-sample VCF files. 

[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/yaacoo/multiPGS_py)

## Motivation

Current available tools that calculate PGS still require many manual steps (e.g., to account for flips), require a very strict format, or only accept a cohort VCF, instead of individual genomes (single-sample indexed VCF files). 
This method can score up to 5 different PGS on an individual genome and can easily be applied to many VCF files using a simple bash script or as part of a workflow/job scheduler. 

## Low-memory for efficient parallelization

The program reads both the PGS file and the VCF file one line at a time, and only keeps a PGS dictionary in memory, which is no more than 60-70Mb. This allows parallelization over many VCF files without running into memory issues. 


## What we account for

1. Strand flips (when the DNA strand orientation differs between the reference and the sample)
2. Allele flips (including beta flip to -beta)
3. Strand + allele flips and ambiguous variants
4. After these adjustments, the calculation is simple:
```math
PGS_{individual} = \sum_{i=1}^{n} (dosage_i \times \beta_i) 
```

## What you still need to check and verify

1. Filter/QC of the imputation quality, filter by max(genotype probability) if needed.
2. That you know the VCF genome build, and which column in the PGS file corresponds to the relevant genome build. This should be in the PGS header, and you should edit line 10 of the Python file if needed.
3. That your VCF is bgzipped and indexed by tabix, having the index file in the same path and prefix (e.g. S001.vcf.gz + S001.vcf.gz.tbi)

## How to use

```
python multiPGS_py.py <single sample VCF file> <PGS catalog file/s, up to 5 files>
```

Example:
```
python multiPGS_py.py sample.vcf.gz PGS000001.txt.gz PGS000002.txt.gz
```

* Don't forget to check line 7 of the Python file to verify the correct column name for your genome build!
* Make sure the VCF has a dosage (DS) field.
* The output file is a single line text file (for each VCF) and can easily be concatenated.
* Tested against Plink2 and pgsc_calc, and provided similar results.

## Parallel processing of multiple genomes simultaneously (experimental)

Set the number of CPU cores available per VM (default n_cpu=20) and process genomes in batches based on the number of available CPU cores. 
After each batch is processed, the next batch is loaded. <br> 
Using a cluster of 10 VMs, each with 20 CPU cores and 1 GB of RAM per CPU, you can expect to process 1,000 genomes in less than 45 minutes. <br>

```
python parallel_pgs.py <vcf_list.txt> <pgs_file1.txt.gz> [<pgs_file2.txt.gz> ...]
```


## Dependencies

python 3.6.8 <br>
pysam 0.16.0.1 <br>
numpy 1.18.5

## Limitations of PGS to be aware of

1. PGS are not always meaningful for an individual, but their rank/percentile/Z-score across a population can indicate risk groups.
2. PGS are not always transferrable between different population ancestries. Please read the paper that published the PGS to PGS Catalog to understand its limitations.

 
