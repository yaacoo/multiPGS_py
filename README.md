# multiPGS_py: Fast, simple and low-memory PGS scoring

multiPGS_py is a fast, simple, and low-memory Python method to calculate polygenic scores (PGS/PRS) from any PGS Catalog weights (pgscatalog.org) and single-sample VCF files. 

## Motivation

Current available tools that calculate PGS still require many manual steps (e.g., to account for flips), require a very strict format, or only accept a cohort VCF, instead of individual genomes (single-sample indexed VCF). 
This method can score up to 5 different PGS on an individual genome and can easily be applied to many VCF files using a simple bash script or as part of a workflow/job scheduler. 

## Low-memory for efficient parallelization

The program reads both the PGS file and the VCF file one line at a time, and only keeps a PGS dictionary in memory, which is no more than 60-70Mb. This allows parallelization over many VCF files without running into memory issues. 
Need to make 

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
2. That you know the VCF genome build, and which column in the PGS file corresponds to the relevant genome build. This should be in the PGS header, and you should edit line 7 of the Python file if needed.
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

## Dependencies

Python 3.6.8
pysam 0.16.0.1
numpy 1.18.5

## Limitations of PGS to be aware of

1. PGS are not always meaningful for an individual, but their rank/percentile/Z-score across a population can indicate risk groups.
2. PGS are not always transferrable between different population ancestries. Please read the paper that published the PGS to PGS Catalog to understand its limitations.

 
