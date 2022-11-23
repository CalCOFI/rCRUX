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

# Statement of Need

A major challenge for metabarcoding is the lack of comprehensive, locus-specific reference databases needed for taxonomic assignment. Although curated databases for a few widely used loci exist (e.g. @quast2012silva), most loci lack dedicated reference databases, especially as the number of target loci grows. Thus, there is a need for an effective tool to generate comprehensive databases. Current approaches to generating reference databases are problematic. For example, approaches that rely on in silico (e.g. @boyer2016obitools) PCR fail to capture reference sequences containing primer sequences while approaches that rely on keyword searches are sensitive to inaccurate metadata (e.g. @bengtsson2018metaxa2). To address this some promising solutions like CRUX (Creating reference libraries using existing tools, @Curd:2019) from the Anacapa Toolkit can generate comprehensive reference databases. However, CRUX software relies on a suite of dependencies making it cumbersome to install and is highly computationally intensive making it difficult to run without access to high performance computing. Combined, these issues highlight the need for a user-friendly, comprehensive reference database generating tool to enhance metabarcoding taxonomic assignment efforts.

# Acknowledgements
Dedicated to the late, great Jesse Gomer. Coding extraordinaire and dear friend.

This work benefited from the amazing input of many including Lenore Pipes, Sarah Stinson, Gaurav Kandlikar, and Maura Palacios Mejia.

Support for the development of this tool was provided by CalCOFI, NOAA, Landmark College, and VBRN.

# References
