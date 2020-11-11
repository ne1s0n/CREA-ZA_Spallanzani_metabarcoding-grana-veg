library(plyr)
load_one_OTU_file = function(infile){
  cnames = c(
    'reads', 'domain', 'clade_phylum', 'class', 
    'order', 'family', 'subfamily', 'tribe', 
    'subtribe', 'genus', 'species')
  df = read.table(infile, sep = '\t', header = FALSE, stringsAsFactors = FALSE, col.names = cnames, fill = TRUE, na.strings = '')
  
  #taking care of anomalous unassigned rows
  aur = subset(df, is.na(domain))
  if (nrow(aur) > 0){
    #removing the anomalous line
    df = subset(df, !is.na(domain))
    #adding back to the unassigned bin
    df[df$domain == 'Unassigned', 'reads'] = df[df$domain == 'Unassigned', 'reads'] + sum(aur$reads)
  }
  return(df)
}

load_all_OTU_files = function(infolder){
  res = list()
  for (f in list.files(path = infolder, pattern = 'derep.nochim.centroids_table.krona.tsv', full.names = TRUE, recursive = TRUE)){
    #recovering the sample name from filename, we have a couple of options
    sample_name = basename(f)
    sample_name = gsub(sample_name, pattern = '-Rub.erne_maxee*', replacement = '')
    sample_name = gsub(sample_name, pattern = '.erne_maxee.*', replacement = '')
    sample_name = gsub(sample_name, pattern = '.derep.nochim.*', replacement = '')
    
    #loading the data
    res[[sample_name]] = load_one_OTU_file(f)
  }
  
  return(res)
}

stats_table = function(raw_OTUs_list){
  res = NULL
  for (sample in names(raw_OTUs_list))
    res = rbind(res, data.frame(
      sample = sample,
      reads = sum(raw_OTUs_list[[sample]]$reads),
      unassigned = subset(raw_OTUs_list[[sample]], domain == 'Unassigned')$reads,
      alpha = nrow(raw_OTUs_list[[sample]]) - 1
    ))
  
  warning('add info on region/province/type here')
  
  res$unassigned_ratio = res$unassigned / res$reads
  return(res)
}

build_OTU_table = function(raw_OTU_list, level='clade_phylum'){
  res = data.frame(row.names = names(raw_OTU_list))
  
  for(sample in names(raw_OTU_list)){
    #compacting for the required level
    curr = raw_OTU_list[[sample]]
    curr = curr[!is.na(curr[,level]),]
    curr = ddply(curr, level, function(x){
      return(data.frame(reads = sum(x$reads)))
    })
    
    #adding the results
    for (i in nrow(curr)){
      res[sample, curr[i, level]] = curr[i, 'reads']
    }
  }
  
  return(res)
}

load_sample_codes = function(){
  samples = read.csv('~/research/CREA-ZA_Spallanzani_metabarcoding-grana-veg/data/private/Transcodifica Campioni_cleaned.csv', stringsAsFactors = FALSE, row.names = 1)
  samples = samples[,c("Id_Lab", "type", "caseificio", "MESE", "data_prelievo")]
  samples$province = sapply(samples$caseificio, FUN = substr, start=1, stop=2)
  samples[samples$province == '', 'province'] = '-'
  
  return(samples)
}
infolder = '~/research/CREA-ZA_Spallanzani_metabarcoding-grana-veg/data/private/delivery_20200706/'
delivery_20200706 = load_all_OTU_files(infolder)
infolder = '~/research/CREA-ZA_Spallanzani_metabarcoding-grana-veg/data/private/delivery_20201009/'
delivery_20201009 = load_all_OTU_files(infolder)
tot = c(delivery_20200706, delivery_20201009)

OTU = build_OTU_table(tot, 'clade_phylum') 
OTU = build_OTU_table(tot, 'class') 
OTU = build_OTU_table(tot, 'order') 
OTU = build_OTU_table(tot, 'family') 
OTU = build_OTU_table(tot, 'genus') 
OTU = build_OTU_table(tot, 'species') 


# 'reads', 'domain', 'clade_phylum', 'class', 
# 'order', 'family', 'subfamily', 'tribe', 
# 'subtribe', 'genus', 'species')

