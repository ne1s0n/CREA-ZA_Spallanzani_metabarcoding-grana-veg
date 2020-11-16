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
  
  #Unassigned, NA and similar all need to be glued in a conventional empty symbol
  df[is.na(df)] = 'UN'
  df[df == 'Unassigned'] = 'UN'
  
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

build_OTU_table = function(raw_OTU_list, level='clade_phylum', min_reads = 5){
  OTU = data.frame(row.names = names(raw_OTU_list))
  
  #---OTU table, absolute values
  for(sample in names(raw_OTU_list)){
    #compacting for the required level
    curr = raw_OTU_list[[sample]]
    curr = curr[!is.na(curr[,level]),]
    curr = ddply(curr, level, function(x){
      return(data.frame(reads = sum(x$reads)))
    })
    
    #if the required level does not pass the filter on reads it gets binned
    #to unassigned
    if(any(curr$reads < min_reads)){
      UN_sel = curr$domain == 'UN'
      bad_sel = curr$reads < min_reads
      curr[UN_sel, 'reads'] = curr[UN_sel, 'reads'] + curr[bad_sel, 'reads']
      curr = curr[!bad_sel,]
    }
    curr = subset(curr, reads >= min_reads)
    
    
    #adding the results
    for (i in 1:nrow(curr)){
      OTU[sample, curr[i, level]] = curr[i, 'reads']
    }
  }
  
  #---relative OTUs count
  #tot reads per sample
  rps = rowSums(OTU, na.rm = TRUE)
  
  OTU_rel = OTU
  for (i in 1:nrow(OTU)){
    OTU_rel [i,] = OTU [i,] / rps[i]
  }
  
  #---relative OTUs count, but for screen printing
  OTU_rel_screen = round(OTU_rel * 100, digits = 3)
  
  #---metadata
  meta = load_sample_codes()
  meta = meta[rownames(OTU),]
  
  return(list(
    OTU = OTU,
    OTU_rel = OTU_rel,
    OTU_rel_screen = OTU_rel_screen,
    meta = meta
  ))
}

load_sample_codes = function(){
  samples = read.csv('~/research/CREA-ZA_Spallanzani_metabarcoding-grana-veg/data/private/Transcodifica Campioni_cleaned.csv', stringsAsFactors = FALSE, row.names = 1)
  samples = samples[,c("Id_Lab", "Label_IGA", "type", "caseificio", "MESE", "data_prelievo")]
  samples$province = sapply(samples$caseificio, FUN = substr, start=1, stop=2)
  samples[samples$province == '', 'province'] = '-'
  samples$Label_IGA = gsub(samples$Label_IGA, pattern = '_', replacement = '-')
  rownames(samples) = samples$Label_IGA
  return(samples)
}

filter_OTU = function(OTU, sample_selector){
  #filtering out the unwanted stuff
  OTU$meta    = OTU$meta[sample_selector,]
  OTU$OTU     = OTU$OTU[sample_selector,]
  OTU$OTU_rel = OTU$OTU_rel[sample_selector,]
  OTU$OTU_rel_screen = OTU$OTU_rel_screen[sample_selector,]
  
  #removing empty cols
  empty_col = colSums(is.na(OTU$OTU)) == nrow(OTU$OTU)
  OTU$OTU     = OTU$OTU[,!empty_col, drop=FALSE]
  OTU$OTU_rel = OTU$OTU_rel[,!empty_col, drop=FALSE]
  OTU$OTU_rel_screen = OTU$OTU_rel_screen[,!empty_col, drop=FALSE]
  
  #done
  return(OTU)
}


