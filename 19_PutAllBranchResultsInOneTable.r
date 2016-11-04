l## Author: Eric Fontanillas
## Last modification: 30122011
## Subject: Extract stats for Free Ratio Branch model run (distri w inter locus for each branch + median + stdev // plot NdN vs SdS

DIR = '18_statsOutput_perBranches_TABLES/'

l<-list.files(path=DIR)

len <- length(l)

tab <- numeric()

for (i in 1:len){
  branch = l[i]
  #print(branch)
  file_IN <- paste(DIR,branch,sep="")
  print(file_IN)
  dat <- read.table(file_IN,  header = TRUE, sep=";")

  row.names(dat) <- branch
  
  #dat2 <- cbind(branch,dat)
  #print(dat)
  tab <- rbind(tab,dat)
  
}
print(tab)

write.csv(tab, "20_table.csv")
