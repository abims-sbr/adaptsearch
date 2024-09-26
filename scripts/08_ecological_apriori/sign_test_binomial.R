#!/usr/bin/env Rscript

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

for(file in list_files){
    name_OUT = paste(path_OUT,"/", file, sep='')
    name_IN = paste(path_IN,"/", file, sep='')

    a <- read.table(name_IN,sep=",", header=T,row.names=1)

    # Pp vs HOTS
    success_Pp <- a$Pp_vs_HOTS[[1]]
    trial_Pp <- a$Pp_vs_HOTS[[4]]
    b_Pp<- binom.test(success_Pp, trial_Pp)
    pv_Pp <- b_Pp$p.value[[1]]  ## p.value
    estimate_Pp <- b_Pp$estimate[[1]] ## probability of success
    conf_low_Pp <- b_Pp$conf.int[[1]] ## conf interval (lower value)
    conf_up_Pp <- b_Pp$conf.int[[2]] ## conf interval (upper value)

    # Pg vs HOTS
    success_Pg <- a$Pg_vs_HOTS[[1]]
    trial_Pg <- a$Pg_vs_HOTS[[4]]
    b_Pg<- binom.test(success_Pg, trial_Pg)
    pv_Pg <- b_Pg$p.value[[1]]  ## p.value
    estimate_Pg <- b_Pg$estimate[[1]] ## probability of success
    conf_low_Pg <- b_Pg$conf.int[[1]] ## conf interval (lower value)
    conf_up_Pg <- b_Pg$conf.int[[2]] ## conf interval (upper value)
    
    # Ap vs COLDS
    success_Ap <- a$Ap_vs_COLDS[[1]]
    trial_Ap <- a$Ap_vs_COLDS[[4]]
    b_Ap<- binom.test(success_Ap, trial_Ap)
    pv_Ap <- b_Ap$p.value[[1]]  ## p.value
    estimate_Ap <- b_Ap$estimate[[1]] ## probability of success
    conf_low_Ap <- b_Ap$conf.int[[1]] ## conf interval (lower value)
    conf_up_Ap <- b_Ap$conf.int[[2]] ## conf interval (upper value)
    
    # Ps vs COLDS
    success_Ps <- a$Ps_vs_COLDS[[1]]
    trial_Ps <- a$Ps_vs_COLDS[[4]]
    b_Ps<- binom.test(success_Ps, trial_Ps)
    pv_Ps <- b_Ps$p.value[[1]]  ## p.value
    estimate_Ps <- b_Ps$estimate[[1]] ## probability of success
    conf_low_Ps <- b_Ps$conf.int[[1]] ## conf interval (lower value)
    conf_up_Ps <- b_Ps$conf.int[[2]] ## conf interval (upper value)

    # join LOWER and UPPER values for confidence interval
    #CI_Pp <- paste(conf_low_Pp, conf_up_Pp)
    #CI_Pg <- paste(conf_low_Pg, conf_up_Pg)
    #CI_Ap <- paste(conf_low_Ap, conf_up_Ap)
    #CI_Ps <- paste(conf_low_Ps, conf_up_Ps)
    
    # Create rows to add to the table "a"
    p_value <- c(pv_Pp,pv_Pg,pv_Ap,pv_Ps)
    estimate <- c(estimate_Pp,estimate_Pg,estimate_Ap,estimate_Ps)
    confidence_interval_lower <- c(conf_low_Pp,conf_low_Pg,conf_low_Ap,conf_low_Ps)
    confidence_interval_upper <- c(conf_up_Pp,conf_up_Pg,conf_up_Ap,conf_up_Ps)
    
    a <- rbind(a,p_value)
    a <- rbind(a,estimate)
    a <- rbind(a,confidence_interval_lower)
    a <- rbind(a,confidence_interval_upper)

    write.table(a,file=name_OUT,sep=",")
}

