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
    affiliation: "1, 2"
  - name: Zack Gold
    orcid: 0000-0003-0490-7630
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
 - name: California Cooperative Oceanic Fisheries Investigations, Scripps Institution of Oceanography & NOAA Southwest Fisheries Science Center, USA
   index: 2
 - name: Southern California Coastal Watershed Research Project
   index: 3
 - name: Universidad Autónoma de Madrid, Spain
   index: 4
 - name: Vermont Biomedical Research Network, University of Vermont, USA
   index: 5

date: 22 November 2022
bibliography: paper.bib

Zenodo-doi: 10.5281/zenodo.7349493
---

# Summary

eDNA metabarcoding is increasingly used to survey biological communities using common universal and novel genetic loci. There is a need for an easy to implement computational tool that can generate metabarcoding reference libraries for any locus, and are specific and comprehensive. We have reimagined CRUX (@Curd:2019) and developed the rCRUX package for the R system for statistical computing (@R_language_2021) to fit this need by generating taxonomy and fasta files for any user defined primer set. The typical workflow involves using *in silico* PCR (e.g. @ye2012primer) to acquire a set of matching sequences containing metabarcode primer sequences. These sequences or "seeds" recovered from the *in silico* PCR step are then used as references used in a second database query to capture complementary sequence that lack one or both primers. This search step, is used to iteratively align seed sequences against a local NCBI database for matches. To increase speed while maintaining comprehensiveness, we implement a taxonomic rank based stratified random sampling approach. The blast_seeds() step results in a comprehensive database of primer specific reference barcode sequences from NCBI. Lastly, the derep_and_clean_db() is used to de-replicate the database by DNA sequence whereby identical sequences are collapsed into a representative sequence read. If there are multiple possible taxonomic paths for a sequence read, the taxonomic path is collapsed to the lowest taxonomic agreement. The fasta and taxonomy output format is accepted by and can easily be formatted for most taxonomic classifiers.

# Acknowledgements
Dedicated to the late, great Jesse Gomer. Coding extraordinaire and dear friend.

This work benefited from the amazing input of many including Lenore Pipes, Sarah Stinson, Gaurav Kandlikar, and Maura Palacios Mejia.

Support for the development of this tool was provided by CalCOFI, NOAA, Landmark College, and VBRN.

# References
