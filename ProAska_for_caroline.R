
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

dms <- remove.poor.quality.snps(d2, min_repro=0.96, max_missing=0.3) %>% sample.one.snp.per.locus.random(., seed=214241)

dms_askania <- remove.by.list(dms, dms$sample_names[!is.na(dms$meta$analyses[,"Analysis.1"])]) %>%
  remove.poor.quality.snps(., min_repro=0.96, max_missing=0.3) 


#### Stats ####
#can be used to remove specific samples

stats <- species_site_stats(dms = dms, pop_var = "sp", missing = 0.3, maf = 0.05,
                            site_var = "site", remove_monomorphic = TRUE)
# write.xlsx(stats, '/Users/eilishmcmaster/Documents/ProAska/ProAska/outputs/Pro_sp_site_stats.xlsx')
# stats_sp <- multispecies_stats(dms = dms, var = dms$meta$analyses[,"sp"], missing = 0.3, maf = 0.05)
# write.xlsx(stats_sp, '/Users/eilishmcmaster/Documents/ProAska/ProAska/outputs/Pro_sp_stats.xlsx')

#### private alleles ####

fitz_allele_list <- make_allele_list(dms_askania, dms_askania$meta$analyses[,"Analysis.1"], min_af = 0.05)

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

######## kinship #############


kin <- individual_kinship_by_pop(dms_askania, RandRbase, species, dataset,
                                 dms_askania$meta$analyses[,"bouddi_vs_all"],
                                 maf=0.1, mis=0.2, as_bigmat=TRUE)

kin2 <- as.data.frame(kin)

col_fun2 = colorRamp2(c(0,0.25,0.45), c("white", "red","black"))

kin2$sample <- rownames(kin2)

hm_sites2 <- merge(kin2, m2[,c("sample","site","sp", "lat", "long")],
                   by="sample", all.x=TRUE, all.y=FALSE)
hm_sites2 <- hm_sites2[match(rownames(kin2),hm_sites2$sample),]
rownames(hm_sites2) <- hm_sites2[,"sample"]

hm_sites2[,"sample"] <- NULL

#create annotations
aska_site_colours <-  named_list_maker(dms_askania$meta$site, "Set3", 11)

site_ann <- HeatmapAnnotation(Site = hm_sites2$site,
                              col=list(Site=aska_site_colours))

site_ann_right <- rowAnnotation(Site = hm_sites2$site,
                              col=list(Site=aska_site_colours),
                              annotation_name_gp = gpar(fontsize = 0))

hma <- Heatmap( as.matrix(hm_sites2[ , c(1:(nrow(hm_sites2)))]), 
                col=col_fun2, 
                bottom_annotation=c(site_ann),
                right_annotation = c(site_ann_right),
                name = "PLINK kinship", #title of legend
                row_names_gp = gpar(fontsize = 6),
                column_names_gp = gpar(fontsize = 6),
                row_names_max_width = unit(15, "cm"),
                border_gp = gpar(col = "black", lty = 1),
                # column_order=order(hm_sites2$lat),
                # row_order=order(hm_sites2$lat)
)


draw(hma, merge_legend = TRUE)

# Set the file name and parameters
filename <- "/Users/eilishmcmaster/Documents/ProAska/ProAska/outputs/plots/kinship.png"
width <- 16
height <- 12
dpi <- 600
units <- "cm"

# Set up the PNG device
png(filename, width = width, height = height, units = units, res = dpi)

# Draw the plot
draw(hma, merge_legend = TRUE)

# Turn off the PNG device
dev.off()


#### relatives ####

# Finding the clones
#https://kateto.net/netscix2016.html


kin_pl2 <- round(kin,2)
kin_pl2[lower.tri(kin_pl2, diag=TRUE)] <- NA
kin_long_pl <- melt(kin_pl2)
kin_long_pl <- kin_long_pl[kin_long_pl$value>0 & !is.na(kin_long_pl$value),]
kin_long_pl <- kin_long_pl[kin_long_pl$Var1!=kin_long_pl$Var2,]
colnames(kin_long_pl)[3] <-"kin"
