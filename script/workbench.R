setwd('~/research/CREA-ZA_Spallanzani_metabarcoding-grana-veg/')
source('script/data_manager.R')

level = 'genus'
level = 'species'

# DATA LOAD ---------------------------------------------------------------

infolder = '~/research/CREA-ZA_Spallanzani_metabarcoding-grana-veg/data/private/delivery_20200706/'
delivery_20200706 = load_all_OTU_files(infolder)
infolder = '~/research/CREA-ZA_Spallanzani_metabarcoding-grana-veg/data/private/delivery_20201009/'
delivery_20201009 = load_all_OTU_files(infolder)
raw_OTUs = c(delivery_20200706, delivery_20201009)

#building the OTU table
OTU_milk = build_OTU_table(raw_OTUs, level) 

#subset to non-control and remove bad sample
sel = OTU_milk$meta$type != 'ctrl' 
OTU_milk = filter_OTU(OTU_milk, sel)

#  ALPHA DIVERSITY --------------------------------------------------------
#View(OTU_milk$meta)
library(vegan)
library(ggplot2)
library(reshape2)
alpha = data.frame(
  Shannon = diversity(OTU_milk$OTU_noUN, index = 'shannon'),
  Simpson = diversity(OTU_milk$OTU_noUN, index = 'simpson'),
  #invsimpson = diversity(OTU_milk$OTU_noUN, index = 'invsimpson'),
  type = OTU_milk$meta$type, 
  sample = rownames(OTU_milk$meta)
)

#moving to long format, for ggplot
alpha_long = melt(alpha, id.vars = c('sample', 'type'), variable.name = 'alpha')

p = ggplot(alpha_long, aes(x=alpha, y=value, fill=type)) + geom_boxplot() + 
  theme(legend.position = 'bottom', legend.title = element_blank(), axis.title = element_blank()) +
  ggtitle(paste('Alpha diversity at rank:', level))
p




# SHANNON DIVERSITY -------------------------------------------------------
library(vegan)
#View(OTU_milk$OTU_noUN)
OTU_milk$meta$shannon_div = diversity(OTU_milk$OTU_noUN)
sha = diversity(OTU_milk$OTU_noUN)

library(ggplot2)
p = ggplot(OTU_milk$meta, aes(x=type, y=shannon_div)) + geom_boxplot()
p


#shapiro test for normality
grana_shannon = subset(OTU_milk$meta, type == 'grana')$shannon_div
nongrana_shannon = subset(OTU_milk$meta, type == 'non-grana')$shannon_div

grana_shap = shapiro.test(grana_shannon)
nongrana_shap = shapiro.test(nongrana_shannon)

hist(grana_shannon)
hist(nongrana_shannon)

res.ftest <- var.test(shannon_div ~ type, data = OTU_milk$meta)
res.ftest
t.test(x=grana_shannon, y=nongrana_shannon)


#ergo: grana diversity is non normal. I cannot use t test
#I'll use non parametric two-samples Wilcoxon rank test.

w = wilcox.test(grana_shannon, nongrana_shannon)

#not statistically different, but look at the histograms!

df = subset(OTU_milk$meta, type == 'grana')
model = lm(shannon_div ~ province, data = df)
anova(model)
summary(model)


# BETA DIVERSITY ----------------------------------------------------------
library(vegan)
library(ggplot2)

#single level
level = 'genus'
level = 'family'
level = 'species'

#building the OTU table
OTU_milk = build_OTU_table(raw_OTUs, level) 

#subset to non-control and remove bad sample
sel = OTU_milk$meta$type != 'ctrl' 
OTU_milk = filter_OTU(OTU_milk, sel)




beta_dist = vegdist(OTU_milk$OTU_noUN, index = "bray")  
mds <- metaMDS(beta_dist, try = 100, trymax = 100)
mds_data <- as.data.frame(mds$points)

mds_data$code <- rownames(mds_data)
mds_data <- dplyr::left_join(mds_data, OTU_milk$meta)
mds_data$ccf = grepl(mds_data$Id_Lab, pattern = 'ccf')


mds_data$month = substr(mds_data$MESE, start = 1, stop = 3)
mds_data$month = factor(mds_data$month, levels = c('', 'GEN', 'FEB', 'MAR', 'APR', 'MAG', 'GIU', 'LUG', 'AGO', 'SET', 'OTT', 'NOV', 'DIC'))

shapes_prov = c(
  BS=17,
  CN=18,
  CR=19,
  PC=0,
  TN=1,
  VI=2,
  VR=3
)

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = type)) +
  geom_point()
ggplot(subset(mds_data, province!= '-'), aes(x = MDS1, y = MDS2, color = province, shape=province)) +
  scale_shape_manual(values = shapes_prov) +
  geom_point()
ggplot(subset(mds_data), aes(x = MDS1, y = MDS2, color = ccf, shape=ccf)) +
  geom_point()

#multi levels
levels = c('family', 'genus', 'species')

beta = NULL

for (level in levels){
  #building the OTU table
  OTU_milk = build_OTU_table(raw_OTUs, level) 
  
  #subset to non-control and remove bad sample
  sel = OTU_milk$meta$type != 'ctrl' 
  OTU_milk = filter_OTU(OTU_milk, sel)
  
  beta_dist = vegdist(OTU_milk$OTU_noUN,
                      index = "bray")  
  
}



















beta_dist = vegdist(OTU_milk$OTU_noUN, index = "bray")  
mds <- metaMDS(beta_dist)
mds_data <- as.data.frame(mds$points)

mds_data$code <- rownames(mds_data)
mds_data <- dplyr::left_join(mds_data, OTU_milk$meta)


library(ggplot2)
ggplot(mds_data, aes(x = MDS1, y = MDS2, color = type)) +
  geom_point()
ggplot(subset(mds_data, province!= '-'), aes(x = MDS1, y = MDS2, color = province, shape=province)) +
  geom_point()
