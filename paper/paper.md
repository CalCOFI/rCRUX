---
title: 'rCRUX: Generate CRUX metabarcoding reference libraries in R'
tags:
  - R
  - environmental DNA
  - primer blastn
  - metabarcode
  - in silico PCR
  - taxonomy
  - reference libraries
authors:
  - name: Luna Gal
    affiliation: 1
  - name: Zack Gold
    orcid: 0000-0000-0000-0000
    affiliation: "2, 3"
  - name: Ramon Gallego
    affiliation: 4
  - name: Emily Curd
    orcid: 0000-0003-0336-6852
    affiliation: "1, 5"
    corresponding: true
affiliations:
 - name: Landmark College, USA
   index: 1
 - name: California Cooperative Oceanic Fisheries Investigations, USA
   index: 2
 - name: Pacific Marine Environmental Laboratory, USA
   index: 3
 - name: Universidad Autónoma de Madrid, Spain
   index: 4
 - name: Vermont Biomedical Research Network, University of Vermont, USA
    index: 5

date: 16 November 2022
bibliography: paper.bib

Zenodo-doi: 10.3847/xxxxx
---

# Summary

eDNA metabarcoding is increasingly used to survey biological communities using common universal and novel genetic loci. There is a need for an easy to implement computational tool that can generate metabarcoding reference libraries for any locus, and are specific and comprehensive. We have reimagined CRUX (@Curd:2019) and developed the rCRUX package for the R system for statistical computing (r2018r) to fit this need by generating taxonomy and fasta files for any user defined locus. The typical workflow involves using get_seeds_local() or get_seeds_remote() to simulate in silico PCR to acquire a set of sequences analogous to PCR products containing metabarcode primer sequences. The sequences or "seeds" recovered from the in silico PCR step are used to search databases for complementary sequence that lack one or both primers. This search step, blast_seeds() is used to iteratively align seed sequences against a local NCBI database for matches using a taxonomic rank based stratified random sampling approach. This step results in a comprehensive database of primer specific reference barcode sequences from NCBI. Using derep_and_clean_db(), the database is de-replicated by DNA sequence where identical sequences are collapsed into a representative read. If there are multiple possible taxonomic paths for a read, the taxonomic path is collapsed to the lowest taxonomic agreement.


# Acknowledgements
Dedicated to the late, great Jesse Gomer. Coding extraordinaire and dear friend.

This work benefited from the amazing input of many including Lenore Pipes, Sarah Stinson, Gaurav Kandlikar, and Maura Palacios Mejia.

Support for the development of this tool was provided by CalCOFI, NOAA, Landmark College, and VBRN.

# References
