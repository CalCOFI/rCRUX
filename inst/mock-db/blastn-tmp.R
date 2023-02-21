tmp_fasta <- tempfile()

writeLines(
  c('>A1234 Vibrio harveyi',
    'gcctaacacatgcaagtcgagcggaaacgagttatctgaaccttcggggaacgataacggcgtcgagcggcggacgggtgagtaatgcctaggaaattgccctgatgtgggggataaccattggaaacgatggctaataccgcataatacctwcgggtcaaagagggggaccttcgggcctctcgcgtcaggatatgcctaggtgggattagctagttggtgaggtaatggctcaccaaggcgacgatccctagctggtctgagaggatgatcagccacactggaactgagacacggtccagactcctacgggaggcagc'), 
  tmp_fasta)

readLines(tmp_fasta)

system2('blastn',
        args = 
          c('-db inst/mock-db/blastdb/mock-db',
            '-query ', tmp_fasta,
            '-outfmt "7 delim=, sseqid staxid"'
          )
)