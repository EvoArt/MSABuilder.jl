# Reference alignments from MUSCLE 3.8.31 (Bioconductor msa package; its
# default clustering is "upgma", so MUSCLE's own default "upgmb" is set), for
# each input *.fa: <name>_it1.fa (maxiters = 1), <name>_it2.fa (maxiters = 2)
# and <name>_it16.fa (defaults). Run from the package root after
# test/fixtures/make_inputs.jl.
suppressMessages(library(msa))
dir <- "test/fixtures"
inputs <- list.files(dir, pattern = "^[a-z]+[.]fa$", full.names = TRUE)
for (f in inputs) {
  s <- Biostrings::readDNAStringSet(f)
  for (it in c(1, 2, 16)) {
    m <- msaMuscle(s, maxiters = it, order = "input", type = "dna", cluster = "upgmb")
    out <- sub("[.]fa$", sprintf("_it%d.fa", it), f)
    Biostrings::writeXStringSet(as(unmasked(m), "DNAStringSet"), out, width = 20000)
  }
}
cat(packageDescription("msa")$Version, "\n")
