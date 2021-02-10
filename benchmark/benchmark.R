library(ape)
library(seqinr)
library(tidysq)
library(Biostrings)
library(bioseq)
library(dplyr)
library(pbapply)

set.seed(981099)

### reading fasta format benchmark

generate_dna_ex <- function(n, len, alph) 
  tidysq::write_fasta(random_sq(n, len, "dna_bsc", 0.1*len), 
              paste0(1:n, alph), 
              paste0("benchmark/dna_ex_n", n, "_l", len, "_a", length(alph), ".fasta"))

alphs <- list("dna_bsc")
ns <- c(10, 100, 1000, 10000)
lens <- c(10, 100, 1000, 10000)

lapply(ns, function(n) {
  lapply(lens, function(len) {
    lapply(alphs, function(alph) {
      generate_dna_ex(n, len, alph)
    })
  })
})

f_read <- list(tidysq = function(x) tidysq::read_fasta(x, "dna_bsc", safe_mode = FALSE),
               seqinr = function(x) seqinr::read.fasta(x), 
               ape = function(x) ape::read.FASTA(x), 
               Biostrings = function(x) Biostrings::readDNAStringSet(x),
               bioseq = function(x) bioseq::read_fasta(x, "DNA"))

f_char <- list(tidysq = function(x) as.character(x[["sq"]]),
               seqinr = function(x) seqinr::getSequence(x), 
               ape = function(x) as.character(x), 
               Biostrings = function(x) sapply(x, toString),
               bioseq = function(x) as.character(x))

f_cons <- list(tidysq = function(x) tidysq::sq(x, "dna_bsc", safe_mode = FALSE),
               seqinr = function(x) lapply(x, function(s) seqinr::as.SeqFastadna(seqinr::s2c(s))),
               ape = function(x) ape::as.DNAbin(x), 
               Biostrings = function(x) Biostrings::DNAStringSet(x),
               bioseq = function(x) bioseq::new_dna(x))

f_tran <- list(tidysq = function(x) tidysq::translate(x[["sq"]]),
               seqinr = function(x) seqinr::getTrans(x),
               ape = function(x) ape::trans(x),
               Biostrings = function(x) Biostrings::translate(x),
               bioseq = function(x) bioseq::seq_translate(x))

num_pkgs <- 5
num_funcs <- 4

results <- do.call(rbind, pblapply(1:10, function(dummy) {
  do.call(rbind, lapply(ns, function(n) {
    do.call(rbind, lapply(lens, function(len) {
      do.call(rbind, lapply(alphs, function(alph) {
        do.call(rbind, lapply(names(f_read), function(i) {
          elapsed_time_r <- system.time(seq_from_fasta <- f_read[[i]](paste0("benchmark/dna_ex_n", n, "_l", len, "_a", 
                                                                             length(alph),".fasta")))
          elapsed_time_char <- system.time(seq_string <- f_char[[i]](seq_from_fasta))
          elapsed_time_cons <- system.time(seq_from_string <- f_cons[[i]](seq_string))
          elapsed_time_tran <- system.time(seq_trans <- f_tran[[i]](seq_from_fasta))
          
          data.frame(package = i,
                     alph_size = length(alph), 
                     sq_len = len,
                     num_sq = n, 
                     type = c("read", "char", "cons", "tran"),
                     file_size = file.size(c(paste0("dna_ex_n", n, "_l", len, "_a", 
                                                    length(alph), ".fasta"))),
                     obj_size = c(as.numeric(object.size(seq_from_fasta)),
                                  as.numeric(object.size(seq_string)),
                                  as.numeric(object.size(seq_from_string)),
                                  as.numeric(object.size(seq_trans))),
                     time_value = c(unname(elapsed_time_r[3]),
                                    unname(elapsed_time_char[3]),
                                    unname(elapsed_time_cons[3]),
                                    unname(elapsed_time_tran[3]))
          )
        }))
      }))
    }))
  }))
}))

alph <- alphs[[1]]

results$file_size <- file.size(rep(
  unlist(
    lapply(
      sapply(ns, function(n) paste0("benchmark/dna_ex_n", n)),
      function(t) paste0(t, "_l", lens, "_a", length(alph), ".fasta"))),
  each = num_funcs * num_pkgs))



write.csv(results, "benchmark/results.csv", row.names = FALSE)
