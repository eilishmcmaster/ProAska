
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
library(vegan)
library(data.table)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

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


#### FST ####


# FST and distance calculations #####################################################################
# remove samples that are NA in species column
fst_dms <- remove.by.list(dms, dms$sample_names[!is.na(dms$meta$analyses[,"species"])]) 

sppop_freq <- as.data.frame(table(fst_dms$meta$analyses[,'species']))
not_n1_speciess <- as.vector(sppop_freq[sppop_freq$Freq<=5,1]) #remove groups where n<=5
not_n1_samples <- fst_dms$sample_names[which(!(fst_dms$meta$analyses[,'species'] %in% not_n1_speciess))]
fst_dms <- remove.by.list(fst_dms, not_n1_samples)


# calculate FST and geodist
gds_file <- dart2gds(fst_dms, RandRbase, species, dataset)
pFst      <- population.pw.Fst(fst_dms, fst_dms$meta$analyses[,'species'], RandRbase,species,dataset, maf_val=0.05, miss_val=0.2) #calculates genetic distance 
pS        <- population.pw.spatial.dist(fst_dms, fst_dms$meta$analyses[,'species']) #calculates geographic distance between populations


####plot IBD plot

# Make self comparisons NA
diag(pFst$Fst) <- NA
diag(pS$S) <- NA

#Mantel test 
man <- mantel(xdis = pS$S, ydis = pFst$Fst, permutations = 10000, na.rm = TRUE) #mantel test, finds if matrices are signficantly similar
man

# mantel plot
Fst_sig <- cbind(melt(pS$S), unlist(as.list(pFst$Fst)))
colnames(Fst_sig)[3] <- "Geo_dist"
colnames(Fst_sig)[4] <- "Fst"
Fst_sig$Geo_dist2 <-Fst_sig$Geo_dist/1000 

meta_agg <- m2 %>%
  group_by(sp, species) %>%
  summarize(lat = mean(lat, na.rm=TRUE),
            long = mean(long,na.rm=TRUE),
            .groups = 'drop')

# adding metadata for speciess
Fst_sig2 <- merge(Fst_sig, distinct(meta_agg[,c("species","sp")]), by.x="Var1", by.y="species", all.y=FALSE)
Fst_sig2 <- merge(Fst_sig2, distinct(meta_agg[,c("species","sp")]), by.x="Var2", by.y="species", all.y=FALSE)
Fst_sig2$same_sp <- ifelse(Fst_sig2$sp.x == Fst_sig2$sp.y, "Within group", "Between group")

fstp1<- ggplot(Fst_sig2, aes(x= Geo_dist2, y=Fst, color=same_sp))+geom_point(size=1, alpha=0.3)+
  labs(x="Distance (km)", y="FST", colour="Comparison")+
  # facet_zoom(x=Geo_dist2<2, zoom.size=1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="bottom")
fstp1

ggsave("ProAska/outputs/plots/ProAska_manning_fst.png",
       fstp1, width = 15, height = 15, units = "cm", dpi=600)

paste("Mantel statistic r is", round(man$statistic, 3), ", P =", man$signif)


# Make heatmaps
# geo dist
geo_d <-pS$S #this is a square matrix
mat_dist <- geo_d/1000 # convert to km 

#FST
mat_fst <-pFst$Fst
agg <- unique(m2[, c("species", "sp")]) # create aggregated df of species and species
mat_dist <- merge(mat_dist, agg, by.x=0, by.y="species", all.y=FALSE) #add aggregated df to mat_fst (fst)
rownames(mat_dist) <- mat_dist$Row.names

mat_dist$Row.names <- NULL
mat_dist <- mat_dist[match(colnames(mat_dist)[1:nrow(mat_dist)],rownames(mat_dist)),]

order_hm <- Heatmap(mat_fst,
                    cluster_rows = TRUE,
                    cluster_columns = TRUE)
od <- colnames(mat_fst)[column_order(order_hm)]

mat_fst = mat_fst[od, od]
mat_dist = mat_dist[od, c(od,"sp")]

# specify fst heatmap colours 
gene_col <-  colorRamp2(c(0,0.5,1), c("#8DD3C7", "white", "#FB8072"))


#specify geo heatmap colours
palette <-  colorRamp2(c(0, max(mat_dist[,1:nrow(mat_dist)], na.rm=TRUE)), c("white", "#80B1D3"))

library(RColorBrewer)
species_colours <-  named_list_maker(dms$meta$analyses[,'sp'], "Spectral", 11)

row_Subpopulation_ann <- rowAnnotation(Species = mat_dist$sp,
                                       col=list(Species=species_colours),
                                       na_col="white",
                                       annotation_legend_param = list(labels_gp=gpar(fontsize=6),#fontface="italic",
                                                                      title_gp=gpar(fontsize=8)),
                                       annotation_name_gp = gpar(fontsize = 0),
                                       annotation_name_side="top")

bottom_Subpopulation_ann <- HeatmapAnnotation(Species = mat_dist$sp, 
                                              col = list(Species = species_colours),
                                              annotation_name_gp = gpar(fontsize = 0),
                                              show_legend = FALSE,
                                              annotation_name_side="right",
                                              na_col = "white")

geo <- Heatmap(as.matrix(mat_dist[,1:nrow(mat_dist)]),rect_gp = gpar(type = "none"),
               width = nrow(mat_dist)*unit(4, "mm"),
               height = nrow(mat_dist)*unit(4, "mm"),
               col=palette,na_col="white",
               bottom_annotation = bottom_Subpopulation_ann,
               column_names_gp = gpar(fontsize = 6),
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               name="Distance (km)",
               heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                                           labels_gp = gpar(fontsize = 6)),
               # cluster_rows = TRUE, 
               # cluster_columns = TRUE,
               cell_fun = function(j, i, x, y, w, h, fill) {
                 if(i >= j) {
                   grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                 }
               }
)

# make fst heatmap
gene <- Heatmap(as.matrix(mat_fst[,1:nrow(mat_fst)]), rect_gp = gpar(type = "none"),
                width = nrow(mat_fst)*unit(4, "mm"),
                height = nrow(mat_fst)*unit(4, "mm"),
                right_annotation = row_Subpopulation_ann,
                col=gene_col,na_col="grey",
                row_names_gp = gpar(fontsize = 6),
                column_names_gp = gpar(fontsize = 0),
                border_gp = gpar(col = "black", lty = 1),
                name="FST",
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                                            labels_gp = gpar(fontsize = 6)),
                cell_fun = function(j, i, x, y, w, h, fill) {
                  if(i <= j) {
                    grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                    grid.text(sprintf("%.2f", mat_fst[,1:nrow(mat_fst)][i, j]), x, y, gp = gpar(fontsize = 4))
                  }
                })

gene_width <- nrow(mat_fst)*unit(4, "mm")

draw(geo + gene, ht_gap = -gene_width)

# Set the file name and parameters
filename <- "/Users/eilishmcmaster/Documents/ProAska/ProAska/outputs/plots/fst_plot_overlay2.png"
width <- 16
height <- 12
dpi <- 600
units <- "cm"

# Set up the PNG device
png(filename, width = width, height = height, units = units, res = dpi)

# Draw the plot
draw(geo + gene, ht_gap = -gene_width)

# Turn off the PNG device
dev.off()



#### map ####
#filter meta here 
meta_map <- m2
# get map bounds
bound <- c(
  left = min(meta_map$long, na.rm = T) - 0.05, bottom = min(meta_map$lat, na.rm = T) - 0.05,
  right = max(meta_map$long, na.rm = T) + 0.05, top = max(meta_map$lat, na.rm = T) + 0.05)

# get stadia map
ggmap::register_stadiamaps("api")

map <- ggmap::get_stadiamap(bbox = bound, zoom=8, scale=4, 
                            maptype="stamen_terrain",
                            color="color") %>%
  ggmap::ggmap()


map+coord_fixed()+
  geom_point(aes(x = long, y = lat, color = sp),
             data = meta_map, size = 2, alpha = 1)+theme_few()

# map + coord_fixed()+
#   geom_point(aes(x = long, y = lat, color = pop_large, shape = pop_large),
#              data = meta_all_fitz, size = 2, alpha = 1) + 
#   scale_colour_manual(values=pop_colours, na.translate=TRUE,
#                       guide = guide_legend("Subpopulation"))+ #override.aes = list(shape = 21)
#   scale_shape_manual(values=pop_large_shapes, na.translate=FALSE,
#                      guide = guide_legend("Subpopulation"))+    theme(legend.position = "bottom") +
#   theme(legend.title=element_blank()) +
#   labs(x = "Longitude", y = "Latitude") +
#   ggsn::scalebar(meta_all_fitz, dist = 1, dist_unit = "km", location = "bottomright", 
#                  st.bottom = F, st.size = 3, st.dist = 0.04,border.size =0.5,
#                  transform = TRUE, model = "WGS84", height = 0.02) +
#   annotation_north_arrow(location = "tr", which_north = "true", pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"), style = north_arrow_fancy_orienteering)+
#   xlim(min(meta_all_fitz$long)-0.005, max(meta_all_fitz$long)+0.005)+
#   ylim(min(meta_all_fitz$lat)-0.005, max(meta_all_fitz$lat)+0.005)+
#   theme_bw()
