library(tidyverse)
library(tidysq)
library(ggseqlogo)
library(PepTools)

pep_dat <- read_tsv(file = "https://raw.githubusercontent.com/leonjessen/tensorflow_rstudio_example/master/data/ran_peps_netMHCpan40_predicted_A0201_reduced_cleaned_balanced.tsv")
pep_dat_sq <- mutate(pep_dat, peptide = construct_sq(peptide, type = "ami"))

pep_dat %>% filter(label_chr=='SB') %>% pull(peptide) %>% ggseqlogo()
pep_dat_sq %>% filter(label_chr=='SB') %>% pull(peptide) %>% as.character %>% ggseqlogo()

pep_dat %>% filter(label_chr=='SB') %>% head(1) %>% pull(peptide) %>% pep_plot_images
pep_dat_sq %>% filter(label_chr=='SB') %>% head(1) %>% pull(peptide) %>% as.character %>% pep_plot_images

pep_dat %>% filter(data_type == 'train') %>% pull(peptide) %>% pep_encode
pep_dat %>% filter(data_type == 'train') %>% pull(peptide) %>% as.character %>% pep_encode

