
#####Load required packages#####
library(RRtools)
library(dartR)
library(cowplot)#
library(geosphere)#
library(tiff)#
library(maptools)#
library(oz)#
library(mapplots)#
library(ggplot2)#
library(ggrepel)#
library(poppr)#
library(SNPRelate)
library(gdsfmt)
library(openxlsx)
library(readxl)
library(ggthemes)
library(data.table)

# load(paste0("C:/Users/cristc01/Documents/ProAska/PSFsaved_scripts4autom.R")) #update 
devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/main/sos_functions.R?raw=TRUE")


maindir = ""
RandRbase <- ""#"C:/Users/cristc01/Documents/"

species <- "ProAska"
dataset <- "DPro23-8557" #this is the identifier provided by DArT. it's in the filenames they provide
m2 <- custom.read(species, dataset) #read custom metadata csv

# import data
d1 <- new.read.dart.xls.onerow(RandRbase,species,dataset,topskip, nmetavar, euchits=FALSE, altcount = TRUE) #read DArT file
# import meta
mm1 <- read.meta.data.full.analyses.df(d1, RandRbase, species, dataset)
# attach meta
d2 <- dart.meta.data.merge(d1, mm1)
# keep only pro aska samples
d3 <- remove.by.list(d2, d2$sample_names[!is.na(d2$meta$analyses[,"Analysis.1"])]) #can be used to remove specific samples 

dx <- remove.poor.quality.snps(d3, min_repro=0.96, max_missing=0.2) %>% sample.one.snp.per.locus.random(., seed=214241)


stats <- species_site_stats(dms = dx,pop_var = "sp",missing = 0.2, maf = 0.05,site_var = "Analysis.1",remove_monomorphic = TRUE)

fitz_allele_list <- make_allele_list(dx, dx$meta$analyses[,"Analysis.1"], min_af = 0.05)

private_total_alleles <- calculate_private_alleles(fitz_allele_list)
private_total_alleles

stats_and_pa <- merge(stats, private_total_alleles, by.x="site",by.y="population", all.y=TRUE)

library(ggVennDiagram)

# compare alleles of gv and sw versus all other sites
fitz_allele_list2 <- list("All other\nsites"=(unlist(fitz_allele_list[c(1,3,4,5,6,8)],use.names = FALSE))%>% as.vector(),
                          "Glenning\nValley"=fitz_allele_list[2]%>% unlist() %>% as.vector(),
                          "Strickland\nWest"=fitz_allele_list[7]%>% unlist()%>% as.vector())

# fitz_allele_list2 <- list("All other\nsites"=(unlist(fitz_allele_list[c(1,3,4,5,6,2,8)],use.names = FALSE))%>% as.vector(),
#                           # "Glenning\nValley"=fitz_allele_list[2]%>% unlist() %>% as.vector(),
#                           "Strickland\nWest"=fitz_allele_list[7]%>% unlist()%>% as.vector())

# fitz_allele_list3 <- list("all"=unique(unlist(fitz_allele_list[c(1,3,5,6,8)],use.names = FALSE))%>% as.vector(),
#                           "gv"=fitz_allele_list[2]%>% unlist(),
#                           "sw"=fitz_allele_list[7]%>% unlist(),
#                           "np"=fitz_allele_list[4]%>% unlist())

ggVennDiagram(fitz_allele_list, label_alpha = 0, edge_size = 0.5, label_size=4)

ggVennDiagram(fitz_allele_list2, label_alpha = 0, edge_size = 0.5, label_size=4)+
  scale_fill_gradient(low="lightblue",high = "white", trans="log")+theme(legend.position = "none")


# ggVennDiagram(fitz_allele_list3, label_alpha = 0, edge_size = 0.5, label_size=4)+
#   scale_fill_gradient(low="lightblue",high = "white", trans="log")+theme(legend.position = "none")

