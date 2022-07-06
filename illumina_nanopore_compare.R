
##### 
##Relative abundance of Illumina dataset - no rarefaction
setwd("~/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/16S_Nanopore/Seaquencing_all_fastq_pass")
asv = read.csv("../../Illumina_reads/run_separation/seaquencing_ASVseq_taxa_4analysis.csv")

library(reshape2)
library(tidyverse)

asv = asv %>% filter(Kingdom == "Bacteria")


abund = asv[,2:41]
sums = colSums(abund)
abund2 = t(t(abund)/sums*100)
abund3 = data.frame(Phylum = asv$Phylum, abund2) 
abund4 = abund3 %>% group_by(Phylum) %>% summarize_all(list(sum))
abund4$Phylum[is.na(abund4$Phylum)] <- "Unknown"

#for calculating average abundance of "Unknown" classifications
y = abund4 %>% filter(Phylum == "Unknown") 
mean(as.matrix(y[,2:ncol(y)]))
sd(as.matrix(y[,2:ncol(y)]))

tax_df = abund4
colnames(tax_df) = gsub("JZO.2021.Condor.", "", colnames(tax_df))
colnames(tax_df) = gsub(".Univ.20220222", "", colnames(tax_df))
colnames(tax_df) = gsub("TSU.2021.Condor.", "", colnames(tax_df))
colnames(tax_df) = gsub(".20211213", "", colnames(tax_df))

colnames(tax_df)[2:ncol(tax_df)] = paste("illumina",  colnames(tax_df)[2:ncol(tax_df)], sep = "-")

#write.csv(tax_df, "../../Illumina_reads/run_separation/Illumina_Phylum_summary.csv")




###########################################


#nanopore data
nano = read.csv("Seaquences-no-rarefaction/Phylum_summary.csv",header = TRUE)
colnames(nano) = gsub("barcode[0-9][0-9]_", "", colnames(nano))
colnames(nano) = gsub("_combined", "", colnames(nano))

colnames(nano)[2:ncol(nano)] = paste("nano",  colnames(nano)[2:ncol(nano)], sep = "-")
nano$Phylum[is.na(nano$Phylum)] <- "Unknown"

#join illumina and nanopore tables
tax_comb = full_join(tax_df, nano, by = "Phylum")
tax_comb2 = tax_comb %>% pivot_longer(!Phylum, names_to = "Sample", values_to = "Abundance") %>% separate(Sample, into = c("Tech", "Sample"), sep = "-") 


#change underscores to periods in sample names 
tax_comb2$Sample = gsub("_", ".", tax_comb2$Sample)
tax_comb2$Sample = gsub("2A2.D52", "2A2", tax_comb2$Sample)
tax_comb2$Sample = gsub("2AT.E46", "2AT", tax_comb2$Sample)
tax_comb2$Sample = gsub("2B1.B22", "2B1", tax_comb2$Sample)
tax_comb2$Sample = gsub("2B1.C18", "2B1", tax_comb2$Sample)
tax_comb2$Sample = gsub("51.E21", "51", tax_comb2$Sample)
tax_comb2$Sample = gsub("Transect","", tax_comb2$Sample)
tax_comb2$Sample = gsub("0408", "4.8", tax_comb2$Sample)
tax_comb2$Sample = gsub("04", "0.4", tax_comb2$Sample)
tax_comb2$Sample = gsub("0812", "8.12", tax_comb2$Sample)
tax_comb2$Sample = gsub("1216", "12.16", tax_comb2$Sample)
tax_comb2$Sample = gsub("1620", "16.20", tax_comb2$Sample)
tax_comb2$Sample = gsub("2024", "20.24", tax_comb2$Sample)
tax_comb2$Sample = gsub("2428", "24.28", tax_comb2$Sample)
tax_comb2$Sample = gsub("2832", "28.32", tax_comb2$Sample)
tax_comb2$Sample = gsub("3236", "32.36", tax_comb2$Sample)
tax_comb2$Sample = gsub("3640", "36.40", tax_comb2$Sample)
tax_comb2$Sample = gsub("2430", "24.30", tax_comb2$Sample)
tax_comb2$Sample = gsub("51.Kilo.1.4", "51.Kilo.0.4", tax_comb2$Sample)



tax_comb3 = tax_comb2 %>% pivot_wider(names_from = Tech, values_from = Abundance)
tax_comb3[,3:4][is.na(tax_comb3[,3:4])] <- 0
tax_comb3$Phylum2 = ifelse((tax_comb3$illumina+tax_comb3$nano) > 7, tax_comb3$Phylum, "other")


#colour palette
colours = colorRampPalette(c('brown', 'red',"orange", 'gold',  'forestgreen', 'turquoise', 'lightblue', 'navy', 'purple', 'pink', 'grey', 'black'))(20)

#scatter plot comparison 
gg = ggplot(tax_comb3, aes(x = illumina, y = nano)) + geom_point(aes(colour = Phylum2),size = 2.5) + coord_equal() + geom_smooth(method = "lm") + scale_y_continuous(limits = c(0,85)) + scale_x_continuous(limits = c(0,85)) + labs(x = "Illumina abundance (%)", y = "Nanopore abundance (%)", colour = "Phylum") + scale_colour_manual(values = colours)+theme(legend.key = element_blank(), legend.title = element_text(size = 10), legend.key.height = unit(0.1, 'cm'), panel.border = element_rect(fill = NA, colour = "grey80"),  panel.background = element_blank(), panel.grid.major = element_line(colour = "grey94")) + guides(colour=guide_legend(ncol=1))


#scatter plot comparison log scale 
gg = ggplot(tax_comb3, aes(x = illumina, y = nano)) + geom_point(aes(colour = Phylum2),size = 2.5) + coord_equal() + geom_smooth(method = "lm") + scale_y_continuous(limits = c(NA,85), trans = "log10") + scale_x_continuous(limits = c(NA,85), trans = "log10") + labs(x = "Illumina abundance (%)", y = "Nanopore abundance (%)", colour = "Phylum") + scale_colour_manual(values = colours)+theme(legend.key = element_blank(), legend.title = element_text(size = 10), legend.key.height = unit(0.1, 'cm'), panel.border = element_rect(fill = NA, colour = "grey80"),  panel.background = element_blank(), panel.grid.major = element_line(colour = "grey94")) + guides(colour=guide_legend(ncol=1))

#scatter plot comparison sqrt scale - with 1:1 line
gg = ggplot(tax_comb3, aes(x = illumina, y = nano))+ geom_abline(intercept = 0, slope = 1, colour = "grey30", linetype="dashed") + geom_point(aes(colour = Phylum2),size = 2.5) + coord_equal() + geom_smooth(method = "lm", colour = 'red') + scale_y_continuous(limits = c(NA,85), trans = "sqrt") + scale_x_continuous(limits = c(NA,85), trans = "sqrt") + labs(x = "Illumina abundance (%)", y = "Nanopore abundance (%)", colour = "Phylum") + scale_colour_manual(values = colours)+theme(legend.key = element_blank(), legend.title = element_text(size = 10), legend.key.height = unit(0.1, 'cm'), panel.border = element_rect(fill = NA, colour = "grey80"),  panel.background = element_blank(), panel.grid.major = element_line(colour = "grey94")) + guides(colour=guide_legend(ncol=1))


#summarize by phylum - create ratio column
tax_comb3$Ratio = tax_comb3$illumina/tax_comb3$nano
tax_comb3$Ratio[!is.finite(tax_comb3$Ratio)] <- NA
tax_comb3$Ratio2 = tax_comb3$nano/tax_comb3$illumina
tax_comb3$Ratio2[!is.finite(tax_comb3$Ratio2)] <- NA
tax_comb4 = tax_comb3 %>% group_by(Phylum) %>% summarise(avg_ratio_illum = mean(Ratio, na.rm = TRUE),  avg_illum = mean(illumina), avg_nano = mean(nano), max_illum = max(illumina), max_nano = max(nano))
#tax_comb4$ratio_combo = ifelse(tax_comb4$avg_ratio_illum >1, tax_comb4$avg_ratio_illum, tax_comb4$avg_ratio_nano)
#tax_comb5 <- tax_comb4[order(-tax_comb4$ratio_combo),]\
tax_comb5 <- tax_comb4[order(-tax_comb4$avg_ratio_illum),]
tax_comb5$Phylum <- factor(tax_comb5$Phylum,levels=unique(tax_comb5$Phylum))
tax_comb5 = tax_comb5[complete.cases(tax_comb5[ , 2]),]
tax_comb5 = tax_comb5[tax_comb5$avg_ratio_illum>0,]

xx = ggplot(tax_comb5, aes(x = avg_ratio_illum, y = Phylum))+ geom_vline(xintercept =1, colour = "red", alpha = 0.25) + geom_point(aes(size = avg_illum, colour= avg_nano)) + labs(x = "Illumina abundance : Nanopore abundance", y = "", size = "Illumina \nabundance (%)", colour = "Nanopore \nabundance (%)" ) +theme(axis.text.y = element_text(size = 7), panel.border = element_rect(fill = NA, colour = "grey80"), legend.key = element_blank(), panel.background = element_blank(), panel.grid.major = element_line(colour = "grey99")) + scale_x_continuous(trans = "log10", breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8)) + scale_radius(breaks = c(1,10,20,40), range = c(1,6)) + scale_colour_gradient(high = "#28393E", low = "#90ADB6") 

#max abundance
xx = ggplot(tax_comb5, aes(x = avg_ratio_illum, y = Phylum))+ geom_vline(xintercept =1, colour = "red", alpha = 0.25) + geom_point(aes(size = max_illum, colour= max_nano)) + labs(x = "Illumina abundance : Nanopore abundance", y = "", size = "Max Illumina \nabundance (%)", colour = "Max Nanopore \nabundance (%)" ) +theme(axis.text.y = element_text(size = 7), panel.border = element_rect(fill = NA, colour = "grey80"), legend.key = element_blank(), panel.background = element_blank(), panel.grid.major = element_line(colour = "grey99")) + scale_x_continuous(trans = "log10", breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8)) + scale_radius(breaks = c(1,10,20,40), range = c(1,6)) + scale_colour_gradient(high = "#28393E", low = "#90ADB6")

##box plot with points
tax_comb4 <- tax_comb3[order(-tax_comb3$Ratio),]
tax_comb4$Phylum <- factor(tax_comb4$Phylum,levels=unique(tax_comb4$Phylum))
xx = ggplot(tax_comb4, aes(x = Ratio, y = Phylum))+ geom_vline(xintercept = 1, colour = "red") + geom_boxplot(outlier.shape = NA) + geom_point(alpha = 0.5, size = 0.75) + scale_x_continuous(trans = "log10") + labs(y = "", x = "Illumina abundance:Nanopore abundance") + theme(axis.text.y = element_text(size = 7), panel.border = element_rect(fill = NA, colour = "grey80"), legend.key = element_blank(), panel.background = element_blank(), panel.grid.major = element_line(colour = "grey99"))


#### just choose high abundance phyla
chloro = tax_comb4 %>% filter(Phylum == "Campylobacterota" | Phylum == "Calditrichota" | Phylum == "Proteobacteria" | Phylum == "Caldatribacteriota" | Phylum == "Planctomycetota" | Phylum == "Chloroflexi" | Phylum == "Desulfobacterota" | Phylum == "Bacteroidota" | Phylum == "Unknown" | Phylum == "Latescibacterota" | Phylum == "Acidobacteriota" | Phylum == "Cyanobacteria" | Phylum == "Methylomirabilota" | Phylum == "NB1-j" | Phylum == "Actinobacteriota" | Phylum == "Verrucomicrobiota" | Phylum == "Patescibacteria")

xx = ggplot(chloro, aes(x = Ratio, y = Phylum))+ geom_vline(xintercept = 1, colour = "red") + geom_boxplot(outlier.shape = NA) + geom_point(alpha = 0.5, size = 0.75) + scale_x_continuous(trans = "log10", breaks = c(0.1, 0.25,0.5, 1,2,4,10)) + labs(y = "", x = "Illumina abundance:Nanopore abundance") + theme(axis.text.y = element_text(size = 7), panel.border = element_rect(fill = NA, colour = "grey80"), legend.key = element_blank(), panel.background = element_blank(), panel.grid.major = element_line(colour = "grey99"))


#correlation - pearson - no rarefaction
 cor(tax_comb3$illumina, tax_comb3$nano)
0.9081141

#correlation - spearman - no rarefaction
cor(tax_comb3$illumina, tax_comb3$nano, method = "spearman")
 0.8869522

#no rarefaction
s = lm(formula = tax_comb3$nano ~ tax_comb3$illumina)
summary(s)
Call:
lm(formula = tax_comb3$nano ~ tax_comb3$illumina)

Residuals:
    Min      1Q  Median      3Q     Max 
-21.254  -0.224  -0.224  -0.122  38.228 

Coefficients:
                   Estimate Std. Error t value Pr(>|t|)    
(Intercept)        0.223766   0.043116    5.19 2.26e-07 ***
tax_comb3$illumina 0.845601   0.007424  113.90  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 2.193 on 2758 degrees of freedom
Multiple R-squared:  0.8247,	Adjusted R-squared:  0.8246 
F-statistic: 1.297e+04 on 1 and 2758 DF,  p-value: < 2.2e-16



################
###Make an NMDS based on phyla level abundance

library(vegan)
yy = tax_comb3 %>% select(Phylum, Sample, illumina, nano) %>% pivot_longer(!c(Phylum,Sample), names_to = "Tech", values_to = "Abundance")
yy$Name = paste0(yy$Sample, "_", yy$Tech)
yy2 = yy %>% pivot_wider(id_cols = Name, names_from = Phylum, values_from = Abundance)
yym = as.matrix(yy2[,2:ncol(yy2)])
set.seed(123)

nmds = metaMDS(yym, distance = "bray")
#Stress: 0.09994489 

data.scores = as.data.frame(scores(nmds)$sites)
data.scores$Name = yy2$Name
data.scores = data.scores %>% separate(Name, into = c("Sample", "Tech"), sep = "_")
data.scores$Sample2 = data.scores$Sample
data.scores = data.scores %>% separate(Sample2, into = c("Site", "Subsite", "Depth1", "Depth2"), sep = "\\.")

#coloured by tech
xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + geom_point(aes(colour = Tech), size = 3) + theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey20"), legend.key = element_blank()) + scale_colour_manual(values = c("#ED254E", "#023778"))
 
#coloured by site 
xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + geom_point(aes(colour = Site), size = 3) + theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey20"), legend.key = element_blank()) + scale_colour_manual(values = c("#ED254E", "#023778", "#A9E2A2", "#94A2B3", "#F9DC5C"))
 
 #coloured by subsite 
xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + geom_point(aes(colour = Subsite), size = 3) + theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey20"), legend.key = element_blank()) + scale_colour_manual(values = c("#ED254E", "#023778", "#A9E2A2", "#94A2B3", "#F9DC5C"))
  
#subsite and depth
xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + geom_point(aes(colour = Subsite, size = as.numeric(Depth1))) + theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey20"), legend.key = element_blank()) + scale_colour_manual(values = c("#ED254E", "#023778", "#A9E2A2", "#94A2B3", "#F9DC5C")) + scale_radius(range = c(1,6)) + labs(size = "Depth (cm)")
  
  
 
 
###Anosim and Mantel tests - no rarefaction
ano = anosim(yym, data.scores$Tech, distance = "bray", permutations = 9999)
#ANOSIM statistic R: 0.06384
      #Significance: 0.006

ano = anosim(yym, data.scores$Site, distance = "bray", permutations = 9999)
#ANOSIM statistic R: 0.4657 
     # Significance: 1e-04 

ano = anosim(yym, data.scores$Subsite, distance = "bray", permutations = 9999)
#ANOSIM statistic R: 0.3814 
#Significance: 1e-04 


#Mantel tests
dist.abund = vegdist(yym, method = "bray")
dist.temp = dist(data.scores$Depth1, method = "euclidean")
abund_temp = mantel(dist.abund, dist.temp, method = "spearman", permutations = 9999, na.rm = TRUE)

#Mantel statistic r: 0.3333  
     # Significance: 1e-04 

aa = as.vector(dist.abund)
tt = as.vector(dist.temp)
plot(aa,tt)
at = data.frame(aa = aa, tt = tt)
at2 = ggplot(at, aes(x = aa, y = tt)) + geom_point(alpha = 0.5) + theme_bw() + labs(x = "Bray Curtis Dissimilarity between samples", y = "Depth (cm) between samples") 



########
##Indicator Species 

 tax_comb6 = tax_comb3 %>% select(Phylum, Sample, illumina, nano) %>% pivot_longer(cols = c(illumina, nano), names_to = "tech", values_to = "Abundance")

tax_comb6$Sample2 = paste0(tax_comb6$Sample,"_", tax_comb6$tech)
tax_comb7 = tax_comb6 %>% pivot_wider(names_from = Phylum, id_cols= Sample2, values_from = Abundance)
tax_comb8 = tax_comb7 %>% separate(Sample2, into = c("Sample", "Tech"), sep = "_")
abund = tax_comb8[,3:ncol(tax_comb8)]
tech = tax_comb8$Tech
inv = multipatt(abund, tech, func = "r.g", control = how(nperm=9999))
summary(inv)


###Illumina data set rarefaction
setwd("~/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/16S_Nanopore/Seaquencing_all_fastq_pass")
asv = read.csv("../../Illumina_reads/run_separation/seaquencing_ASVseq_taxa_4analysis.csv")

library(reshape2)
library(tidyverse)

asv = asv %>% filter(Kingdom == "Bacteria")


asvm = melt(asv, id = c("ASVID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))


#library(data.table)
library(splitstackshape)

asvm0 = asvm %>% filter(value != 0)

#error repeating on value column didn't work
#asvm2 = expandRows(asvm0, asvm0$value)
#repeating on 9th column (which is value column did for some reason)
asvm2 = expandRows(asvm0, count = 9)


asvm_split = split(asvm2, asvm2$variable)



#to apply to all elements of asvm list
fun = function(i){
    sub = sample(1:nrow(i), 1000, replace=FALSE)
    x = i[sub,]}

asvm_sub = lapply(asvm_split, fun)
#no rarefaction
#asvm_sub = asvm_split

#apply phyla summary and abundance
fun2 = function(i){
    phy = i %>% filter(Kingdom == "Bacteria") %>% group_by(Phylum) %>% summarise(n = n()) %>% mutate(abund = n/(colSums(as.matrix(n)))*100) %>% select(-n)
    }
phylum_sum = lapply(asvm_sub, fun2)    

tax_df = phylum_sum %>% reduce(full_join, by='Phylum')
colnames(tax_df) = c("Phylum", names(phylum_sum))

#shorten column names 
colnames(tax_df) = gsub("JZO.2021.Condor.", "", colnames(tax_df))
colnames(tax_df) = gsub(".Univ.20220222", "", colnames(tax_df))
colnames(tax_df) = gsub("TSU.2021.Condor.", "", colnames(tax_df))
colnames(tax_df) = gsub(".20211213", "", colnames(tax_df))



colnames(tax_df)[2:ncol(tax_df)] = paste("illumina",  colnames(tax_df)[2:ncol(tax_df)], sep = "-")
tax_df$Phylum[is.na(tax_df$Phylum)] <- "Unknown"
