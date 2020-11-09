setwd('~/research/CREA-ZA_Spallanzani_metabarcoding-grana-veg/')

#read the OTU table
OTU = read.csv('data/private/otu_table-2.csv', stringsAsFactors = FALSE, row.names = 1)
OTU$taxonomy = NULL

#getting the top 30 of most present centroids
centroids_cnt = rowSums(OTU)
OTU = OTU[order(centroids_cnt, decreasing = TRUE),]

#cumulative plot
centroids_cum = centroids_cnt
for (i in 1:length(centroids_cnt)){
  centroids_cum[i] = sum(centroids_cnt[1:i])
}
plot(centroids_cum/sum(centroids_cnt))
(centroids_cum/sum(centroids_cnt))[1000]

#how many reads does the top-100 represents?
writeLines('Percentage of reads in the top-30 most represented centroids')
print(sum(centroids_cnt[1:100]) / sum(centroids_cnt))

#reading the fasta with centroids sequences
fasta = readLines('data/private/centroids/seqs_chimeras_filtered_rep_set.fasta')
reduced_fasta = data.frame(row.names = rownames(OTU)[1:100])
for (i in seq(from=2, to=length(fasta), by=2)){
  #centroid name
  cn = gsub(fasta[i-1], pattern = '>', replacement = '')
  cn = strsplit(cn, split = ' ')[[1]][1]
  
  if (cn %in% rownames(OTU)[1:100]){
    reduced_fasta[cn,'label'] = fasta[i-1]
    reduced_fasta[cn,'seq'] = fasta[i]
  }
}

#writing the selected centroids to a new file
outfile = '~/research/CREA-ZA_Spallanzani_metabarcoding-grana-veg/data/private/centroids/seqs_chimeras_filtered_rep_set.TOP100.fasta'
fp = file(outfile, open = 'a')
for (i in 1:nrow(reduced_fasta)){
  writeLines(reduced_fasta[i, 'label'], con=fp)
  writeLines(reduced_fasta[i, 'seq'], con=fp)
}
close(fp)

