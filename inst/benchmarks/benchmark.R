library(ape)
library(seqinr)
library(tidysq)
library(Biostrings)


generate_dna_ex <- function(n, len, alph) {
  name <- paste0(1:n, "dna")

  char_vec <- unlist(lapply(1L:n, function(i) {
    s <- paste0(sample(alph, len, replace = TRUE), collapse = "")
    paste0(">", name[i], "\n", s, "\n")
  }))

  writeLines(text = char_vec, con = paste0("dna_ex_n", n, "_l", len, "_a", length(alph), ".fasta"))
}

alphs <- list(c("C", "T", "A", "G"), LETTERS)
ns <- 10^(2:5)
lens <- 10^(1:5)

invisible(lapply(ns, function(n) {
  lapply(lens, function(len) {
    lapply(alphs, function(alph) {
      generate_dna_ex(n, len, alph)
    })
  })
}))

library(dplyr)

f_read <- list(tidysq = function(x) tidysq::read_fasta(x, type = "unt"),
               seqinr = function(x) seqinr::read.fasta(x), 
               ape = function(x) ape::read.FASTA(x), 
               Biostrings = function(x) Biostrings::readBStringSet(x))

results <- do.call(rbind, lapply(ns, function(n) {
  do.call(rbind, lapply(lens, function(len) {
    do.call(rbind, lapply(alphs, function(alph) {
      do.call(rbind, lapply(1:length(f_read), function(i) {
        f <- f_read[[i]]
        t0 <- Sys.time()
        s <- f(paste0("dna_ex_n", n, "_l", len, "_a", length(alph),".fasta"))
        t1 <- Sys.time()
        data.frame(package = names(f_read)[i], alph_size = length(alph), sq_len = len, num_sq = n, obj_size = as.numeric(object.size(s)), reading_time = t1-t0)
      }))
    }))
  }))
}))

write.csv(results, "./inst/benchmarks/results.csv", row.names = FALSE)

# library(ggplot2)
# library(reshape2)
# 
# results <- read.csv("./inst/benchmarks/results.csv")
# 
# ggplot(dplyr::filter(results, package != "seqinr"), aes(x = as.factor(num_sq), y = obj_size, fill = package)) +
#   geom_col(position = "dodge") +
#   facet_wrap(~ alph_size + sq_len, labeller = label_both, scales = "free_y")
# 
# ggplot(results, aes(x = as.factor(num_sq), y = reading_time, fill = package)) +
#   geom_col(position = "dodge") +
#   facet_grid(alph_size ~ sq_len, labeller = label_both)
# 
# ggplot(results, aes(x = as.factor(num_sq), y = as.factor(sq_len), fill = reading_time)) +
#   geom_tile() +
#   facet_grid(alph_size ~ package)
