#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

library(tidyverse)

gff_file <- args[1]

reference_start = 1
pdat <- ape::read.gff(gff_file) %>% 
  # as_tibble() %>% 
  select(type, start, end, phase) %>% 
  filter(type == 'CDS') %>% 
  mutate(gene = str_c('cds_', 1:nrow(.)),
         phase = as.integer(as.character(phase)),
         start = start + phase, 
         offset = case_when((start - reference_start + 1 - 0) %% 3 == 0 ~ 0,
                            (start - reference_start + 1 - 1) %% 3 == 0 ~ 1,
                            (start - reference_start + 1 - 2) %% 3 == 0 ~ -1),
         s = (start - reference_start + 1 - offset) / 3)

dat <- map_dfr(1:(nrow(pdat) - 1), function(i) { 
  if ( pdat[i + 1, 4] == 0 ) { int <- filter(pdat, row_number() == i) } 
  if ( pdat[i + 1, 4] == 1 ) { int <- filter(pdat, row_number() == i) %>% mutate(end = end - 2) }
  if ( pdat[i + 1, 4] == 2 ) { int <- filter(pdat, row_number() == i) %>% mutate(end = end - 1) }
  if ( i + 1 == nrow(pdat) ) { int <- bind_rows(int , filter(pdat, row_number() == i + 1)) }
  return(int)
}) %>% 
  mutate(e = ((end - start + 1) / 3) + s - reference_start)  %>% 
  select(gene, start = s, end = e, offset)

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
           