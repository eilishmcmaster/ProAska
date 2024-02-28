

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
# remove loci with missingness >1%
d4 <- remove.poor.quality.snps(d3, min_repro=0.96, max_missing=0.01)
# remove fixed loci
d5 <- remove.fixed.snps(d4)
d6 <- sample.one.snp.per.locus.random(d5, seed=12345)
# remove loci with low minor allele frequency 
d7 <- remove.by.maf(d6, 0.05)

#################### PRIVATE ALLELES ##########################

dx <- remove.poor.quality.snps(d3, min_repro=0.96, max_missing=0.3)%>% sample.one.snp.per.locus.random(., seed=214241)


stats <- species_site_stats(dms = dx,pop_var = "sp",missing = 0.3,maf = 0.05,site_var = "Analysis.1",remove_monomorphic = TRUE)

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

fitz_allele_list3 <- list("all"=unique(unlist(fitz_allele_list[c(1,3,5,6,8)],use.names = FALSE))%>% as.vector(),
                          "gv"=fitz_allele_list[2]%>% unlist(),
                          "sw"=fitz_allele_list[7]%>% unlist(),
                          "np"=fitz_allele_list[4]%>% unlist())

ggVennDiagram(fitz_allele_list, label_alpha = 0, edge_size = 0.5, label_size=4)

ggVennDiagram(fitz_allele_list2, label_alpha = 0, edge_size = 0.5, label_size=4)+
  scale_fill_gradient(low="lightblue",high = "white", trans="log")+theme(legend.position = "none")


ggVennDiagram(fitz_allele_list3, label_alpha = 0, edge_size = 0.5, label_size=4)+
  scale_fill_gradient(low="lightblue",high = "white", trans="log")+theme(legend.position = "none")


##################################  OPTGENMIX ############################################

# ###load required code
# source("/Users/eilishmcmaster/Google Drive/My Drive/Rossetto_lab/sos/RRtools ReCeR/scripts/common_allele_countNARM.r")
# source("/Users/eilishmcmaster/Google Drive/My Drive/Rossetto_lab/sos/RRtools ReCeR/scripts/get_minor_allele_frequenciesNARM.r")
# 
# devtools::install_github("jasongbragg/OptGenMix")
# remotes::install_github("jasongbragg/OptGenMix")
# 
# # Replace 'branch_name' with the specific branch or commit you want to install
# url <- "https://github.com/jasongbragg/OptGenMix.git"
# install.packages(url, repos = NULL, type = "source")

optfiles <- list.files('/Users/eilishmcmaster/Documents/OptGenMix-master/R', full.names = TRUE)

for(i in optfiles){
  source(i)  
}


# library(OptGenMix)
library(phangorn)
library(openxlsx)


set.seed(23908845) #this is for randomisation below to reach "optimum"


# numbers of individuals to optimise for
# N_t_vec <- c(10,12,15,17,20) # collection size
num_steps <- 100
maf <- 0.05

run_optimization_set <- function(gdd,N_t, max_wts,min_wts, initial_weights) { # source("run_optimization_set.r")
  # maximize genetic diversity
  d_opt <- optimize_single_objective(sm=as.matrix(gdd), N_t=N_t, initial_weights=initial_weights, 
                                     weights_max=max_wts, measure="matrix_weighted_mean", max_steps=num_steps, 
                                     max_t=0.01, p_depends_delta=TRUE, weights_min = min_wts)
  
  return(list(N_t=N_t, d_opt=d_opt))
}



N_t_vec_list <- list("Analysis.1"=c(10,12,15,17,20),
                     "strickland_west_only"=c(2,3,4,5),
                     "no_bouddi"=c(2,3,4,5))



analyses <- c("Analysis.1", "strickland_west_only", "no_bouddi")

for(subset in analyses){
  ####
  # d8 <- remove.by.list(d3.1,d3.1$sample_names[(d3.1$meta$analyses[,"Analysis.1"])!="askania_StricklandWest"])
  d8 <- remove.by.list(d7, d7$sample_names[!is.na(d7$meta$analyses[,subset])])%>%
    remove.poor.quality.snps(., min_repro=0.96, max_missing=0.9)%>%
    remove.by.maf(., 0.05)
    
  ####
  dms <- d8
  N_t_vec <- N_t_vec_list[subset][[1]]
  # get the genotype data from filtered dms
  gta <- dms$gt
  # exsitu_samples_index <- which(rownames(gta) %in% exsitu_samples)
  gt_sw_comp <- gta
  i_ub <- c(1:length(rownames(gta))); print(i_ub)
  
  max_wts <- rep(1, nrow(gta)) #the maximum each individual can be chosen
  min_wts <- NULL
  # min_wts <- rep(0, nrow(gt_sw_comp)) #the minimum number of times each individual can be included
  # min_wts[exsitu_samples_index] <- 1 # make sure that ex situ samples are force included
  
  sw_out_list <- list()
  
  gdd <- dist(gta) # turn into distance matrix
  
  # make initial weights and run ioptimisation
  for ( i in 1:length(N_t_vec) ) { # for each collection size
    N_t <- N_t_vec[i]
    cat("\n Running ", N_t, " ...\n")
    
    if(is.null(min_wts)){
      initial_weights = rep(0, nrow(gt_sw_comp))
      rand_samp <- sample(i_ub)[1:N_t] #print(rand_samp)
      initial_weights[rand_samp] <- 1 #(print(initial_weights))
    }else{
      iNt_x <- N_t - length(exsitu_samples_index)
      print(iNt_x)
      initial_weights <- rep(0, nrow(gt_sw_comp))
      # Choose additional individuals other than ex situ samples
      non_exsitu_samples <- setdiff(i_ub, exsitu_samples_index)
      selected_samples <- sample(non_exsitu_samples, size = iNt_x)
      # Set ex situ samples to 1 and selected individuals to 1 in initial_weights
      initial_weights[exsitu_samples_index] <- 1
      initial_weights[selected_samples] <- 1
    }
    sw_out_list[[ i ]] <- run_optimization_set(gdd=gdd,N_t= N_t, max_wts=max_wts, min_wts=min_wts, initial_weights=initial_weights)  
  }
  
  
  # Create a PNG file for saving the plot
  png(paste0(species, "/outputs/ProAska_",subset,"_d_opt_plot.png"), width = 100, height = 100, units="mm", res=300)
  # Generate the plot
  plot(sw_out_list[[1]]$d_opt$value, main = "T_max = 0.01")
  # Close the PNG device to save the file
  dev.off()
  
  out_alleles <- mat.or.vec(length(N_t_vec),2)
  sw_maf    <- get_minor_allele_frequencies(gt_sw_comp) 
  
  threshold_maf <- maf
  i_sw_common <- which(sw_maf > threshold_maf) 
  
  for (i in 1:length(N_t_vec)) {
    sol_vec <- sw_out_list[[i]]$d_opt$weight[num_steps,]
    common_alleles  <- common_allele_count(gt_sw_comp, sol_vec) 
    out_alleles[i,1] <- i
    out_alleles[i,2] <- length(intersect( which(common_alleles[[2]] > 0), i_sw_common))
  }
  
  ##############################################################################
  ### compare to random translocation populations of 
  ### different sizes
  ###note that the random results can be higher than optimum due to 
  rNt <- N_t_vec
  nreps <- num_steps
  rvals <- list()
  sdvals <- list()
  minvals <- list()
  maxvals <- list()
  
  for (i in 1:length(rNt)) {
    ivals <- c()
    iNt <- rNt[i]
    # iNt <- iNt-length(exsitu_samples_index)
    
    for (j in 1:nreps) {
      ran_vec <- rep(0, nrow(gt_sw_comp))
      
      # Choose additional individuals other than ex situ samples
      # non_exsitu_samples <- setdiff(i_ub, exsitu_samples_index)
      selected_samples <- sample(i_ub, size = iNt)
      
      # Set ex situ samples to 1 and selected individuals to 1 in ran_vec
      # ran_vec[exsitu_samples_index] <- 1
      ran_vec[selected_samples] <- 1
      
      common_alleles <- common_allele_count(gt_sw_comp, ran_vec)
      ivals[j] <- length(intersect(which(common_alleles[[2]] > 0), i_sw_common))
    }
    
    mvals <- mean(ivals)
    rvals[i] <- mvals
    sdvals[i] <- sd(ivals)
    minvals[i] <- min(ivals)
    maxvals[i] <- max(ivals)
    
    cat("\n Running ", i, " ...\n")
  }
  
  
  # current_ex_alleles <- common_allele_count(gt_sw_comp[exsitu_samples_index,], w=rep(1, length(exsitu_samples_index))) 
  # 
  # current_ex_alleles2 <- length(intersect( which(current_ex_alleles[[2]] > 0), i_sw_common))/length(i_sw_common)
  
  
  allele_props <- data.frame(n=N_t_vec, 
                             random_mean = (unlist(rvals)),
                             random_sd = unlist(sdvals),
                             random_min = unlist(minvals),
                             random_max = unlist(maxvals),
                             optimised = (out_alleles[,2]),
                             total_common_alleles = length(i_sw_common))
  
  write.xlsx(as.data.frame(allele_props),paste0(species, "/outputs/ProAska_",subset,"_alleleprops.xlsx"))
  
  allele_props_long2 <- reshape2::melt(allele_props,value.name = "common_alleles", 
                                       id.vars=c("n","total_common_alleles","random_sd", "random_min", "random_max"),
                                       variable.name="selection")
  
  allele_props_long3 <- allele_props_long2 %>%
    mutate(random_sd = ifelse(selection == "optimised", NA, random_sd))%>%
    mutate(random_min = ifelse(selection == "optimised", NA, random_min))%>%
    mutate(random_max = ifelse(selection == "optimised", NA, random_max))#%>%
  # mutate(across('selection', str_replace, 'random_mean', 'Random'))%>%
  # mutate(across('selection', str_replace, 'optimised', 'Optimised'))
  
  allele_props_long3$prop_common_alleles <- allele_props_long3$common_alleles/
    allele_props_long3$total_common_alleles
  # allele_props_long3$prop_sd <- as.numeric(allele_props_long3$random_sd)/
  #   allele_props_long3$total_common_alleles
  
  
  allele_prop_plot <- ggplot(allele_props_long3, aes(x=n, y=prop_common_alleles, group=selection, color=selection))+
    theme_few()+
    geom_line(size=0.2, alpha=0.7)+
    # scale_shape_manual(values=c(1,2))+
    labs(y="Proportion of common major alleles", shape="Population", color="Selection")+
    # geom_errorbar(aes(ymin = prop_common_alleles - prop_sd, ymax = prop_common_alleles + prop_sd), width = 0.5, size=0.5)+
    geom_hline(yintercept = 0.9, col="grey", linetype="dotted")+
    geom_hline(yintercept = 0.95, col="grey", linetype="dotted")+
    geom_hline(yintercept = 0.99, col="grey", linetype="dotted")+
    geom_point(size=1)+
    theme(strip.text.y = element_text(angle = 0))+
    theme(legend.position="bottom")
  # geom_hline(yintercept = current_ex_alleles2, col="red")
  
  ggsave(paste0(species, "/outputs/ProtAska_",subset,"_random_vs_opt.png"), allele_prop_plot,
         width=150, height=100, units="mm", dpi=300)
  
  ###get table
  solution_table <- as.data.frame(mat.or.vec(nrow(gt_sw_comp), length(N_t_vec)))
  for (i in 1:length(N_t_vec)) {
    solution_table[,i] <- sw_out_list[[i]]$d_opt$weight[num_steps,] 
  }
  
  
  # solution_table[,1:2] <- c(d4$sample_names, as.character(dms_meta$site))
  colnames(solution_table) <- paste("n=",N_t_vec)
  solution_table$sample <- paste(rownames(gta))
  solution_table <- merge(solution_table, m2[,c("lat", "long", "site","sample")], by="sample",
                          all.y=FALSE)
  
  addrow <- c("Allele proportion", (out_alleles[,2]/length(i_sw_common)), "","","","")
  addrow2 <- c("Random allele proportion", (unlist(rvals)/length(i_sw_common)), "","","","")
  solution_table2 <- rbind(solution_table, addrow, addrow2)
  
  write.xlsx(as.data.frame(solution_table2),paste0(species, "/outputs/ProAska_",subset,"_optgenmix.xlsx"))
  
  
}

  
#   
#              
# ############ missingness all samples #################  
# missingness_gt    <- as.data.frame(rowSums(is.na(dms[["gt"]]))/ncol(dms_rm[["gt"]])) #get missingness of all samples
#   
# write.csv(missingness_gt, "Eido_missingness.csv")
# 
# #######SAM's scripts##########
# ##Load and filter baseline (include all sites) genetic data
# 
# ##Set Max weights – this is where you exclude the individuals you do not want to use
# max_wts <- rep(0, nrow(gta)) #the maximum each individual can be chosen
# i_seed_nh  <- which(dms$meta$analyses[,analysis] != "nursery_seedling")
# max_wts[i_seed_nh]    <- 1
# 
# ###Run Optgenmix
# 
# ###calculate the total common alleles
# #common allel count is 1, minor allele count observed once 
# i_sw_common <- which(sw_maf > threshold_maf) 
# 
# ###get number of common alleles for the optimised selections and calculate portion diversity captured by dividing the numbers by the total 
# for (i in 1:length(N_t_vec)) {
#   sol_vec <- sw_out_list[[i]]$d_opt$weight[num_steps,] #this is the "optimal selection of individuals" 
#   common_alleles  <- common_allele_count(gt_sw_comp, sol_vec) #returns: number_common_alleles=number_common_alleles, minor_allele_counts=minor_allele_counts
#   out_alleles[i,1] <- i
#   out_alleles[i,2] <- length(intersect( which(common_alleles[[2]] > 0), i_sw_common)) #this is the number of common alleles captured by the selection
# }
# 
# 
# 
# 