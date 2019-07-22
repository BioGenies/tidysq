library(ape)
library(seqinr)
library(tidysq)
library(Biostrings)

### reading fasta format benchmark

generate_dna_ex <- function(n, len, alph) {
  name <- paste0(1:n, "dna")

  char_vec <- unlist(lapply(1L:n, function(i) {
    s <- paste0(sample(alph, len, replace = TRUE), collapse = "")
    paste0(">", name[i], "\n", s, "\n")
  }))

  writeLines(text = char_vec, con = paste0("dna_ex_n", n, "_l", len, "_a", length(alph), ".fasta"))
}

alphs <- list(c("C", "T", "A", "G"), LETTERS)
ns <- 10^(1:2)
lens <- 10^(1:3)

invisible(lapply(ns, function(n) {
  lapply(lens, function(len) {
    lapply(alphs, function(alph) {
      generate_dna_ex(n, len, alph)
    })
  })
}))

library(dplyr)

f_read <- list(tidysq = function(x) tidysq::read_fasta(x, type = "unt"),
               ##seqinr = function(x) seqinr::read.fasta(x), 
               ape = function(x) ape::read.FASTA(x), 
               Biostrings = function(x) Biostrings::readBStringSet(x))

set.seed(19296442)
results <- do.call(rbind, lapply(ns, function(n) {
  do.call(rbind, lapply(lens, function(len) {
    do.call(rbind, lapply(alphs, function(alph) {
      do.call(rbind, lapply(1:length(f_read), function(i) {
        do.call(rbind, lapply(1:100, function(dummy) {
          f <- f_read[[i]]
          t0 <- Sys.time()
          s <- f(paste0("inst/benchmarks/dna_ex_n", n, "_l", len, "_a", length(alph),".fasta"))
          t1 <- Sys.time()
          data.frame(package = names(f_read)[i], alph_size = length(alph), sq_len = len, num_sq = n, obj_size = as.numeric(object.size(s)), reading_time = t1-t0)
        }))
      }))
    }))
  }))
}))

res_melted <- results %>%
  select(-obj_size, -alph_size) %>%
  melt(id.vars = c("num_sq", "sq_len", "package")) %>%
  group_by(num_sq, sq_len, package, variable) %>%
  summarise(value = median(value))

ggplot(data = res_melted, aes(x = package, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ sq_len + num_sq, labeller = label_both, scales = "free_y") 

ggplot(data = res_melted, aes(x = sq_len, y = value, color = package)) +
  geom_line() +
  geom_point() +
  facet_grid(variable ~ num_sq)

res_melted_o <- results %>%
  select(-char_to_seq_time, -seq_to_char_time) %>%
  melt(id.vars = c("num_sq", "sq_len", "package")) %>%
  group_by(num_sq, sq_len, package, variable) %>%
  summarise(value = median(value))

ggplot(data = res_melted_o, aes(x = num_sq, y = value, color = package)) +
  geom_line() +
  geom_point() +
  facet_grid(variable ~ sq_len)

ggplot(data = res_melted_o, aes(x = sq_len, y = value, color = package)) +
  geom_line() +
  geom_point() +
  facet_grid(variable ~ num_sq)


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


### converting strings to sequence internal objects

ns <- seq(10, 100, length.out = 4)
lens <- seq(10, 10000, length.out = 10)

f_cons <- list(tidysq = function(x) tidysq::construct_sq(x, type = "unt"),
##            seqinr = function(x) seqinr::as.SeqFastadna(x), 
               ape = function(x) ape::as.DNAbin(strsplit(x, "")), 
               Biostrings = function(x) Biostrings::DNAStringSet(x))

f_char <- list(tidysq = function(x) as.character(x),
##             seqinr = function(x) seqinr::as.SeqFastadna(x), 
               ape = function(x) as.character(x), 
               Biostrings = function(x) sapply(x, toString))

set.seed(19296442)
results <- do.call(rbind, lapply(ns, function(n) {
  do.call(rbind, lapply(lens, function(len) {
    #do.call(rbind, lapply(alphs, function(alph) {
      dna <- sapply(1:n, 
                    function(x) paste0(sample(c("C", "T", "G", "A"), len, replace = TRUE), collapse = ""))
      do.call(rbind, lapply(1:length(f_cons), function(i) {
        do.call(rbind, lapply(1:20, function(dummy) {
          fc <- f_cons[[i]]
          fr <- f_char[[i]]
          t0 <- Sys.time()
          s <- fc(dna)
          t1 <- Sys.time()
          s2 <- fr(s)
          t2 <- Sys.time()
          data.frame(package = names(f_cons)[i], 
                     ###                alph_size = length(alph), 
                     sq_len = len, 
                     num_sq = n, 
                     obj_size = as.numeric(object.size(s)), 
                     char_to_seq_time = t1 - t0,
                     seq_to_char_time = t2 - t1)
        }))
      }))
##  }))
  }))
}))



library(reshape2)

# ggplot(data = results, aes(x = package, y = char_to_seq_time)) +
#   geom_boxplot() +
#   facet_grid(sq_len ~ num_sq, labeller = label_both)
# 
# ggplot(data = results, aes(x = package, y = seq_to_char_time)) +
#   geom_boxplot() +
#   facet_grid(sq_len ~ num_sq, labeller = label_both)

res_melted <- results %>%
  select(-obj_size) %>%
  melt(id.vars = c("num_sq", "sq_len", "package")) %>%
  group_by(num_sq, sq_len, package, variable) %>%
  summarise(value = median(value))

ggplot(data = res_melted, aes(x = package, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ sq_len + num_sq, labeller = label_both, scales = "free_y") 

ggplot(data = res_melted, aes(x = sq_len, y = value, color = package)) +
  geom_line() +
  geom_point() +
  facet_grid(variable ~ num_sq)

res_melted_o <- results %>%
  select(-char_to_seq_time, -seq_to_char_time) %>%
  melt(id.vars = c("num_sq", "sq_len", "package")) %>%
  group_by(num_sq, sq_len, package, variable) %>%
  summarise(value = median(value))

ggplot(data = res_melted_o, aes(x = num_sq, y = value, color = package)) +
  geom_line() +
  geom_point() +
  facet_grid(variable ~ sq_len)

ggplot(data = res_melted_o, aes(x = sq_len, y = value, color = package)) +
  geom_line() +
  geom_point() +
  facet_grid(variable ~ num_sq)
