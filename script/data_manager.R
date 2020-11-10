load_one_OTU_file = function(infile){
  cnames = c(
    'reads', 'domain', 'clade', 'class', 
    'order', 'family', 'subfamily', 'tribe', 
    'subtribe', 'genus', 'species')
  df = read.table(infile, sep = '\t', header = FALSE, stringsAsFactors = FALSE, col.names = cnames, fill = TRUE, na.strings = '')
  return(df)
}

load_all_OTUs = function(infolder){
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


infolder = '~/research/CREA-ZA_Spallanzani_metabarcoding-grana-veg/data/private/delivery_20200706/'
a = load_all_OTUs(infolder)
infolder = '~/research/CREA-ZA_Spallanzani_metabarcoding-grana-veg/data/private/delivery_20201009/'
b = load_all_OTUs(infolder)
