setwd('~/research/CREA-ZA_Spallanzani_metabarcoding-grana-veg/')
infolder = 'data/private/delivery_20201009/raw_reads/'
outfile = 'data/private/delivery_20201009/raw_reads/reads_count.csv'

#for each available infile
res = NULL
all_files = list.files(infolder, pattern = 'fastq.gz$', full.names = TRUE)
for (i in 1:length(all_files)){
  #interface
  writeLines(paste('doing file', i, 'of', length(all_files)))
  
  #extracting line number
  f = all_files[i]
  cmd = paste('zcat', f, '| wc -l')
  lines = as.integer(system(cmd, intern = TRUE))
  
  #storing result
  samplename = gsub(basename(f), pattern = '.fastq.gz$', replacement = '')
  res = rbind(res, data.frame(
    sample = samplename,
    lines = lines,
    reads = lines/4
  ))
}

#saving
write.csv(res, file = outfile, row.names = FALSE)
