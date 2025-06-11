dir = '/Users/dashnowh/Documents/git/inSTRbility'
setwd(dir)
library(ggplot2)
library(cowplot)
# apply cowplot theme
theme_set(theme_cowplot())

# Load the data into a single data frame
# all.data = data.frame()
# for (file in list.files(path = dir, pattern = '*meth.txt', full.names = TRUE)) {
#   this.data = read.table(file, header = F, stringsAsFactors = F)
#   this.sample = gsub('.meth.txt', '', basename(file))
#   this.data$sample = this.sample
#   all.data = rbind(all.data, this.data)
# 
# }

# Understanding methylation values
#meth_bins = data.frame(bin = 0:255, min_prob = (0:255)/256, max_prob = ((0:255) + 1)/256)

metadata = read.csv('FMR1_instability_PureTarget_metadata.tsv', sep = '\t')
#all.data = read.csv('fmr1_instrbility.tsv', sep = '\t')
all.data = read.csv('read-data-allSamples.tsv', sep = '\t')
all.data = merge(all.data, metadata)

all.data = subset(all.data, base_qual > 10)
all.data$allele_length = as.numeric(all.data$allele_length)
all.data$methylated_bases = as.numeric(all.data$methylated_bases)
all.data$umethylated_bases = as.numeric(all.data$umethylated_bases)
all.data$median_meth = as.numeric(all.data$median_meth)
all.data$medianmethlevel = round(all.data$median_meth, 1)
all.data$ismeth = all.data$median_meth >= 0.5
all.data$prop_meth = all.data$methylated_bases/(all.data$methylated_bases + all.data$umethylated_bases)
all.data$prop_meth[is.nan(all.data$prop_meth)] = NA
all.data$prop_meth_level = as.factor(round(all.data$prop_meth, 1))


# for (individual in unique(all.data$Individual)){
#   ggplot(subset(all.data, Individual == individual), aes(x = allele_length, fill = as.factor(medianmethlevel))) + 
#     geom_histogram() + 
#     facet_grid(Region~Individual, scales = 'free') + 
#     scale_fill_manual(values = c('blue4', 'blue3', 'blue2', 'blue1', 'royalblue1', 'snow2', 'tomato', 'red1', 'red2', 'red3', 'red4')) +
#     labs(x = 'Allele length (motifs)', y = 'Reads', fill = 'Median methylation') 
#   ggsave(paste0(individual,'_FMR1_PureTarget_methylation.pdf'))
# }
# ggplot(all.data, aes(x = allele_length, fill = as.factor(medianmethlevel))) + 
#   geom_histogram() + 
#   facet_wrap(Individual~Region, scales = 'free') + 
#   scale_fill_manual(values = c('blue4', 'blue3', 'blue2', 'blue1', 'royalblue1', 'snow2', 'tomato', 'red1', 'red2', 'red3', 'red4')) +
#   labs(x = 'Allele length (motifs)', y = 'Reads', fill = 'Median methylation') 
# ggsave('FMR1_PureTarget_methylation.pdf', width = 20, height = 20)
# 
# ggplot(all.data, aes(x = allele_length, fill = as.factor(medianmethlevel))) + 
#   geom_histogram(aes(y = after_stat(density))) + 
#   facet_grid(Region~Individual, scales = 'free') + 
#   scale_fill_manual(values = c('blue4', 'blue3', 'blue2', 'blue1', 'royalblue1', 'snow2', 'tomato', 'red1', 'red2', 'red3', 'red4')) +
#   labs(x = 'Allele length (motifs)', y = 'Reads', fill = 'Median methylation')  
# ggsave('FMR1_PureTarget_methylation_byregion.pdf', width = 20, height = 20)



ggplot(all.data, aes(x = allele_length, fill = prop_meth_level)) + 
  geom_histogram() + 
  facet_wrap(Individual~Region, scales = 'free') + 
  scale_fill_manual(values = c('blue4', 'blue3', 'blue2', 'blue1', 'royalblue1',
                               'snow2', 'tomato', 'red1', 'red2', 'red3', 'red4'))+
  labs(x = 'Allele length (motifs)', y = 'Reads', fill = 'Prop. Cs in repeat methylated')  
ggsave('FMR1_PureTarget_methylation_prop.pdf', width = 20, height = 20)

for (individual in unique(all.data$Individual)){
  print(individual)
  print(subset(all.data, Individual == individual))
  ggplot(subset(all.data, Individual == individual), aes(x = allele_length, fill = prop_meth_level)) + 
    geom_histogram() + 
    facet_grid(Region~Individual, scales = 'free') + 
    scale_fill_manual(values = c('blue4', 'blue3', 'blue2', 'blue1', 'royalblue1', 'snow2', 'tomato', 'red1', 'red2', 'red3', 'red4')) +
    labs(x = 'Allele length (motifs)', y = 'Reads', fill = 'Prop. Cs in repeat methylated') 
  ggsave(paste0(individual,'_FMR1_PureTarget_methylation_prop.jpg'))
}

ggplot(subset(all.data, allele_length > 400), aes(x = allele_length, fill = prop_meth_level)) +
  geom_histogram() +
  facet_wrap(~Sex) +
  scale_fill_manual(values = c('blue4', 'blue3', 'blue2', 'blue1', 'royalblue1', 'snow2', 'tomato', 'red1', 'red2', 'red3', 'red4')) +
  labs(x = 'Allele length (motifs)', y = 'Reads', fill = 'Prop. Cs in repeat methylated')

ggplot(all.data, aes(x = allele_length, fill = prop_meth_level)) +
  geom_histogram() +
  facet_wrap(~Sex) +
  scale_fill_manual(values = c('blue4', 'blue3', 'blue2', 'blue1', 'royalblue1', 'snow2', 'tomato', 'red1', 'red2', 'red3', 'red4')) +
  labs(x = 'Allele length (motifs)', y = 'Reads', fill = 'Prop. Cs in repeat methylated')

ggplot(subset(all.data, allele_length > 400), aes(x = allele_length, fill = Sex)) +
  geom_histogram() +
  labs(x = 'Allele length (motifs)', y = 'Reads')

ggplot(all.data, aes(x = median_meth, y = prop_meth)) + geom_point()

ggplot(all.data, aes(x = base_qual)) + geom_histogram()
ggplot(all.data, aes(x = base_qual, y = allele_length)) + geom_point() + ylim(0, 100)

ggplot(subset(all.data, Individual == 'FXPM4555'), aes(x = allele_length, fill = prop_meth_level)) + 
  geom_histogram() + 
  facet_grid(Region~Individual, scales = 'free') + 
  scale_fill_manual(values = c('blue4', 'blue3', 'blue2', 'blue1', 'royalblue1', 'snow2', 'tomato', 'red1', 'red2', 'red3', 'red4')) +
  labs(x = 'Allele length (motifs)', y = 'Reads', fill = 'Prop. Cs in repeat methylated') 


#View(subset(all.data, Individual == "CON5497"))
