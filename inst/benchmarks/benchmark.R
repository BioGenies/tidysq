library(ape)
library(seqinr)
library(tidysq)
library(Biostrings)
library(dplyr)
library(pbapply)

set.seed(981099)

### reading fasta format benchmark

generate_dna_ex <- function(n, len, alph) {
  name <- paste0(1:n, "dna")
  
  char_vec <- unlist(lapply(1L:n, function(i) {
    s <- paste0(sample(alph, len + round(rnorm(1, sd = 0.1*len), 0), replace = TRUE), collapse = "")
    paste0(">", name[i], "\n", s, "\n")
  }))
  
  writeLines(text = char_vec, con = paste0("dna_ex_n", n, "_l", len, "_a", length(alph), ".fasta"))
}

alphs <- list(c("C", "T", "A", "G"))
ns <- round(seq(10, 1000, length.out = 5), 0)
lens <- round(seq(10, 1000, length.out = 5), 0)

invisible(lapply(ns, function(n) {
  lapply(lens, function(len) {
    lapply(alphs, function(alph) {
      generate_dna_ex(n, len, alph)
    })
  })
}))

f_read <- list(tidysq = function(x) tidysq::read_fasta(x, type = "unt"),
               seqinr = function(x) seqinr::read.fasta(x), 
               ape = function(x) ape::read.FASTA(x), 
               Biostrings = function(x) Biostrings::readBStringSet(x))

f_cons <- list(tidysq = function(x) tidysq::construct_sq(x, type = "unt"),
               seqinr = function(x) seqinr::as.SeqFastadna(x), 
               ape = function(x) ape::as.DNAbin(x), 
               Biostrings = function(x) Biostrings::DNAStringSet(x))

f_char <- list(tidysq = function(x) as.character(x[["sq"]]),
               seqinr = function(x) seqinr::getSequence(x), 
               ape = function(x) as.character(x), 
               Biostrings = function(x) sapply(x, toString))

results <- do.call(rbind, pblapply(1:20, function(dummy) {
  do.call(rbind, lapply(ns, function(n) {
    do.call(rbind, lapply(lens, function(len) {
      do.call(rbind, lapply(alphs, function(alph) {
        do.call(rbind, lapply(names(f_read), function(i) {
          elapsed_time_r <- system.time(seq_from_fasta <-f_read[[i]](paste0("dna_ex_n", n, "_l", len, "_a", 
                                                                            length(alph),".fasta")))
          elapsed_time_char <- system.time(seq_string <- f_char[[i]](seq_from_fasta))
          elapsed_time_cons <- system.time(seq_from_string <- f_cons[[i]](seq_string))
          
          data.frame(package = i, alph_size = length(alph), 
                     sq_len = len, num_sq = n, 
                     type = c("read", "char", "cons"),
                     obj_size = c(as.numeric(object.size(seq_from_fasta)),
                                  as.numeric(object.size(seq_string)),
                                  as.numeric(object.size(seq_from_string))),
                     time_value_r = c(unname(elapsed_time_r[3]),
                                      unname(elapsed_time_char[3]),
                                      unname(elapsed_time_cons[3]))
          )
        }))
      }))
    }))
  }))
}))

write.csv(results, "./inst/benchmarks/results.csv", row.names = FALSE)
