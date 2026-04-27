# hm-ctDNA 
This repository contains a pipeline for identifying and classifying circulating tumor DNA (ctDNA) fragments based on CpG methylation patterns.
As tumor derived cfDNA is often characterised by distinct hypermethylation patterns and CpG dense regions
We implemented a novel binary classification strategy to distinguish hypermethylated ctDNA (hm-ctDNA) from background cell-free DNA (cfDNA).

## Methodology
Each DNA fragment is evaluated based on:
	Number of CpG sites
	Fraction of methylated CpGs

## Classification Criteria
1. Hypermethylated ctDNA (hm-ctDNA) : A fragment is classified as hypermethylated ctDNA if it contains at least 6 CpG sites, and >75% of CpGs are methylated
2. Non-hypermethylated cfDNA (norm-cfDNA) : A fragment is classified as non-hypermethylated cfDNA if it contains fewer than 6 CpG sites,

## Installation
- >git clone  https://github.com/aishwarya-sekar/hm-ctDNA.git    
- >cd hm-ctDNA

## Usage
python3 ./scripts/run-hm-ctDNa.py

## Requirements
- python3
- samtools
- bedtools
- bismark processed, sorted, indexed, BAM files
- Methyldackel and its dependencies (external [https://github.com/dpryan79/methyldackel])
- ngs-bits and its dependencies (external [https://github.com/imgag/ngs-bits])
