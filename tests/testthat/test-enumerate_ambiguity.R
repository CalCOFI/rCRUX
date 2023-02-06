test_that("enumerate_ambiguity simple works", {
  
  iupac <-
    list( "M" = c("A", "C"),
          "R" = c("A", "G"),
          "W" = c("A", "T"),
          "S" = c("C", "G"),
          "Y" = c("C", "T"),
          "K" = c("G", "T"),
          "V" = c("A", "C", "G"),
          "H" = c("A", "C", "T"),
          "D" = c("A", "G", "T"),
          "B" = c("C", "G", "T"),
          "N" = c("A", "C", "G", "T"),
          "I" = c("A", "T", "C"))
  
  # should generate character vectors of the same values, then we add names
  # and the result should be the same
  enumerated_ambiguties <- lapply(names(iupac), enumerate_ambiguity)
  names(enumerated_ambiguties) <- names(iupac)
  
  expect_identical(enumerated_ambiguties, iupac)
  
})

test_that("enumerate_ambiguity complex works", {
  
  iupac <-
    list( "M" = c("A", "C"),
          "R" = c("A", "G"),
          "W" = c("A", "T"),
          "S" = c("C", "G"),
          "Y" = c("C", "T"),
          "K" = c("G", "T"),
          "V" = c("A", "C", "G"),
          "H" = c("A", "C", "T"),
          "D" = c("A", "G", "T"),
          "B" = c("C", "G", "T"),
          "N" = c("A", "C", "G", "T"),
          "I" = c("A", "T", "C"))
  
  enumerated_ambiguties <- enumerate_ambiguity('MRWS')
  
  # data.frame of combinations
  expected_enumeration_df <- 
    expand.grid(
      iupac$M,
      iupac$R,
      iupac$W,
      iupac$S
    )
  
  # vector of collapsed values, but not in the expected order as enumerated_ambiguties
  expected_enumerations <- apply(expected_enumeration_df, MARGIN = 1, FUN = paste, collapse = '')
  
  # sort on comparison
  expect_equal(sort(enumerated_ambiguties), sort(expected_enumerations))
  
})
