library(ggplot2)
library(reshape2)
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

build_OTU_table = function(raw_OTU_list, level='clade_phylum', min_centroid_reads = 5, min_sample_reads = 1000){
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
    if(any(curr$reads < min_centroid_reads)){
      UN_sel = curr$domain == 'UN'
      bad_sel = curr$reads < min_centroid_reads
      curr[UN_sel, 'reads'] = curr[UN_sel, 'reads'] + curr[bad_sel, 'reads']
      curr = curr[!bad_sel,]
    }
    curr = subset(curr, reads >= min_centroid_reads)
    
    #adding the results
    for (i in 1:nrow(curr)){
      OTU[sample, curr[i, level]] = curr[i, 'reads']
    }
  }
  
  #NAs at this point are considered as zeros
  OTU[is.na(OTU)] = 0
  
  #samples not reaching the minimum amount of reads are removed
  OTU = OTU[rowSums(OTU) >= min_sample_reads,]
  
  #---relative OTUs count
  #tot reads per sample
  rps = rowSums(OTU, na.rm = TRUE)
  
  OTU_rel = OTU
  for (i in 1:nrow(OTU)){
    OTU_rel [i,] = OTU [i,] / rps[i]
  }
  
  #---relative OTUs count, but for screen printing
  OTU_rel_screen = round(OTU_rel * 100, digits = 3)
  
  #---relative OTUs count, ignoring unassigned
  OTU_noUN = OTU
  OTU_noUN$UN = NULL  
  #tot reads per sample
  rps_noUN = rowSums(OTU_noUN, na.rm = TRUE)
  
  OTU_rel_noUN = OTU_noUN
  for (i in 1:nrow(OTU_noUN)){
    OTU_rel_noUN [i,] = OTU_noUN [i,] / rps_noUN[i]
  }
  
  #---relative OTUs count, but for screen printing
  OTU_rel_screen_noUN = round(OTU_rel_noUN * 100, digits = 3)
  
  #---metadata
  meta = load_sample_codes()
  meta = meta[rownames(OTU),]
  meta$caseificio = gsub(meta$caseificio, pattern = '^$', replacement = '[NON-GRANA]')
  meta$code = paste(sep='_', meta$caseificio, meta$Id_Lab, meta$Label_IGA)
  
  #renaming everything
  rownames(meta) = meta$code
  rownames(OTU) = meta$code
  rownames(OTU_rel) = meta$code
  rownames(OTU_rel_screen) = meta$code
  rownames(OTU_noUN) = meta$code
  rownames(OTU_rel_noUN) = meta$code
  rownames(OTU_rel_screen_noUN) = meta$code
  
  #ordering so that production sites are clumped
  neworder = meta$code[order(meta$code)]
  meta = meta[neworder,, drop=FALSE]
  OTU = OTU[neworder,, drop=FALSE]
  OTU_rel = OTU_rel[neworder,, drop=FALSE]
  OTU_rel_screen = OTU_rel_screen[neworder,, drop=FALSE]
  OTU_noUN = OTU_noUN[neworder,, drop=FALSE]
  OTU_rel_noUN = OTU_rel_noUN[neworder,, drop=FALSE]
  OTU_rel_screen_noUN = OTU_rel_screen_noUN[neworder,, drop=FALSE]
  
  return(list(
    OTU = OTU,
    OTU_rel = OTU_rel,
    OTU_rel_screen = OTU_rel_screen,
    meta = meta,
    OTU_noUN = OTU_noUN,
    OTU_rel_noUN = OTU_rel_noUN,
    OTU_rel_screen_noUN = OTU_rel_screen_noUN
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
  OTU$meta    = OTU$meta[sample_selector, , drop=FALSE]
  OTU$OTU     = OTU$OTU[sample_selector, ,drop=FALSE]
  OTU$OTU_rel = OTU$OTU_rel[sample_selector, ,drop=FALSE]
  OTU$OTU_rel_screen = OTU$OTU_rel_screen[sample_selector, , drop=FALSE]
  OTU$OTU_noUN     = OTU$OTU_noUN[sample_selector, , drop=FALSE]
  OTU$OTU_rel_noUN = OTU$OTU_rel_noUN[sample_selector, , drop=FALSE]
  OTU$OTU_rel_screen_noUN = OTU$OTU_rel_screen_noUN[sample_selector, , drop=FALSE]
  
  #removing empty cols
  empty_col = colSums(OTU$OTU) == 0
  OTU$OTU     = OTU$OTU[,!empty_col, drop=FALSE]
  OTU$OTU_rel = OTU$OTU_rel[,!empty_col, drop=FALSE]
  OTU$OTU_rel_screen = OTU$OTU_rel_screen[,!empty_col, drop=FALSE]
  empty_col = colSums(OTU$OTU_noUN) == 0
  OTU$OTU_noUN     = OTU$OTU_noUN[,!empty_col, drop=FALSE]
  OTU$OTU_rel_noUN = OTU$OTU_rel_noUN[,!empty_col, drop=FALSE]
  OTU$OTU_rel_screen_noUN = OTU$OTU_rel_screen_noUN[,!empty_col, drop=FALSE]
  
  #done
  return(OTU)
}

prepare_heatmap = function(mat, min_abundance = 0.01){
  #removing unassigned, if present
  mat$UN = NULL
  
  #removing organisms out of thresholds
  bad = apply(mat, 2, max, na.rm = TRUE) < min_abundance
  mat = mat[, !bad, drop=FALSE]
  
  #sorting by numerosity
  freqsum = colSums(mat)
  mat = mat[, order(freqsum, decreasing = TRUE), drop = FALSE]
  organism_order = colnames(mat)
  
  #building a long format, for plotting
  mat_long = melt(data.frame(mat, sample = rownames(mat)), id.vars = 'sample')
  levels(mat_long$variable) = organism_order
  
  #building a ggplot tile
  gg = ggplot(data=mat_long, aes(y=sample, x=variable, fill=value)) + 
    geom_tile(colour = "gray") +
    scale_fill_gradient(low = "white",high = "blue") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) + 
    #coord_fixed(ratio=1) + 
    theme(
      legend.position = 'bottom',
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=9, face = "italic"),
      axis.text.y = element_text(hjust = 1, size=11, vjust = 0.5),
      axis.title = element_blank()
    )
  
  #returning everything
  return(list(
    data = mat,
    data_long = mat_long,
    gg = gg
  ))
  
}

