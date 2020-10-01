#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

library(tidyverse)

gff_file <- args[1]


gene_start = 1
dat <- ape::read.gff(gff_file) %>%
    #dplyr::as_tibble() %>%
    dplyr::mutate(type = as.character(type),
                  phase = as.numeric(as.character(phase))) %>%
    dplyr::filter(type == 'CDS') %>%
    dplyr::mutate(gene = stringr::str_c('cds_', 1:nrow(.))) %>%
    dplyr::select(gene, start, end, phase) %>% 
    dplyr::mutate(offset_s = if_else((start - gene_start + 1 - 0 + phase) %% 3 == 0, 0, 
                                      if_else((start - gene_start + 1 - 1 + phase) %% 3 == 0, 1, -1)),
                  offset_e = if_else((end - gene_start + 1 - 0) %% 3 == 0, 0, 
                                      if_else((end - gene_start + 1 - 1) %% 3 == 0, 1, 2)),
                  s = (start - gene_start + 1 + offset_s + phase) / 3,
                  e = (end - gene_start + 1 + offset_e) / 3) %>%
    dplyr::select(gene, start = s, end = e, offset = offset_s)

purrr::map(levels(factor(dat$offset)), function(off){
           dat_2 <- dplyr::filter(dat, offset == off)
           if (nrow(dat_2) == 1) {
              dplyr::bind_rows(dat_2, tibble(gene = 'debug',
                                             start = 10000000000,
                                             end = 100000000000,
                                             offset = dat_2$offset)) %>%
              readr::write_csv(paste0('codfreq_gff_offset_', off, '.csv'))
           } else {
              dat_2  %>%
              readr::write_csv(paste0('codfreq_gff_offset_', off, '.csv'))
           }
})
           