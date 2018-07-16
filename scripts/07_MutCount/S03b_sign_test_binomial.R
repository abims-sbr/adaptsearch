#!/usr/bin/env Rscript
#coding: utf-8
#Author : Eric Fontanillas (2010) - Victor Mataigne (2018)

# binom test - null hypothesis : 
# The test computes the number of times where the value of the variable A is higher to the value of the variable B.
# Under the null hypothesis, there is no difference between variables : A > B in 50% of the trials

options(warn = -1)

library(optparse)

option_list = list(
  make_option(c('-i', '--indir'), type='character', default=NULL, help='input directory name', metavar='character'),
  make_option(c('-o', '--outdir'), type='character', default=NULL, help='output directory name', metavar='character')
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
path_IN <- opt$indir
path_OUT <- opt$outdir
files_OUT <- paste(path_OUT,"/*",sep='')
unlink(files_OUT)  ### clean the path_OUT folder from previous files
list_files <- list.files(path_IN)

# iterate over file
for(file in list_files){
    name_OUT = paste(path_OUT,"/", file, sep='')    
    name_IN = paste(path_IN,"/", file, sep='')

    b <- read.table(name_IN, sep=",", header=T, row.names=1)

    supplemental_cells <- matrix(0, nrow=4, ncol=ncol(b))
    colnames(supplemental_cells) <- colnames(b)
    row.names(supplemental_cells) <- c('p-value', 'probability_of_success', 'confidence_interval_low', 'confidence_interval_high')

    # iterate over species groups
    for (i in 1:ncol(b)) {
      if (b[4,i] != 0) {
        binom_test <- binom.test(b[1,i], b[4,i])
        supplemental_cells[1,i] <- binom_test$p.value
        supplemental_cells[2,i] <- binom_test$estimate
        supplemental_cells[3,i] <- binom_test$conf.int[1]
        supplemental_cells[4,i] <- binom_test$conf.int[2]
      } else {
        supplemental_cells[1,i] <- NA
        supplemental_cells[2,i] <- NA      
        supplemental_cells[3,i] <- NA
        supplemental_cells[4,i] <- NA
      }
    }

    final <- rbind(b, supplemental_cells)
    write.csv(final, file=name_OUT)
}

