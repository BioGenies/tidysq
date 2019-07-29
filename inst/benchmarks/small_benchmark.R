library(dplyr)
library(tidysq)
library(pbapply)

alphs <- list(c("C", "T", "A", "G"))
ns <- c(10, 20)
lens <- c(100, 10000)

f_cons <- list(tidysq = function(x) tidysq::construct_sq(x, type = "unt"),
               tidysq2 = function(x) tidysq::construct_sq2(x, type = "unt"),
               tidysq3 = function(x) tidysq::construct_sq3(x, type = "unt"), 
               ape = function(x) ape::as.DNAbin(x), 
               Biostrings = function(x) Biostrings::DNAStringSet(x))

results <- do.call(rbind, pblapply(1:20, function(dummy) {
  do.call(rbind, lapply(ns, function(n) {
    do.call(rbind, lapply(lens, function(len) {
      do.call(rbind, lapply(alphs, function(alph) {
        seq_string <- unlist(lapply(1L:n, function(i) {
          paste0(sample(alph, len + round(rnorm(1, sd = 0.1*len), 0), replace = TRUE), collapse = "")
        }))
        do.call(rbind, lapply(names(f_cons), function(i) {
          elapsed_time_cons <- system.time(seq_from_string <- f_cons[[i]](seq_string))
          
          data.frame(package = i, alph_size = length(alph), 
                     sq_len = len, num_sq = n, 
                     obj_size = c(as.numeric(object.size(seq_from_string))),
                     time_value = c(unname(elapsed_time_cons[3]))
          )
        }))
      }))
    }))
  }))
}))

library(ggplot2)
melt_res <- results %>% 
  group_by(package, sq_len, num_sq) %>% 
  summarise(obj_size = median(obj_size),
            time_value = median(time_value))

filter(melt_res) %>% 
  ggplot(aes(x = sq_len, y = time_value, color = package)) +
  geom_point() +
  geom_line() +
  facet_wrap( ~ num_sq, labeller = label_both, scales = "free_y") +
  ggtitle("cons time") +
  theme_bw() +
  theme(legend.position = "bottom")
