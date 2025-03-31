library(ggplot2)
library(ggpubr)    
library(ggridges)
library(scales) 
library(R.matlab)
library(MatchIt) 
library(xtable)
library(dplyr) 



# modified script based on https://github.com/LenaDorfschmidt/sex_differences_adolescence

# adult cell types
lake <- read.csv("C:\\Users\\lihon\\Desktop\\abcd_study\\pnc\\gene_decoding\\Lake2018.csv")
lake = lake[-which(lake$Cluster=='Ast_Cer'),] 
lake = lake[-which(lake$Cluster=='OPC_Cer'),]
lake = lake[-which(lake$Cluster=='Purk1'),]
lake = lake[-which(lake$Cluster=='Purk2'),]
lake = lake[-which(lake$Cluster=='Per'),]
lake = lake[-which(lake$Cluster=='Gran'),]
lake$Cluster <- as.factor(lake$Cluster)

###
load("C:\\Users\\lihon\\Desktop\\abcd_study\\pnc\\gene_decoding\\sex_differences_adolescence\\Data\\connectivity_matrices\\nspn.main.RData")

###
in.path = "C:\\Users\\lihon\\Desktop\\abcd_study\\pnc\\vis_v2_r400_ro1"
pls1 <- read.csv(paste0(in.path, "\\corr_rbag_gene_pnc.csv")) # Read in corr outputs
colnames(pls1)[1]<- 'hgnc_symbol'

# # Load in available gene lengths and chromosome assignments
pls1.genes = genes %>% left_join(pls1, by='hgnc_symbol')      # Match chromosome assignment and gene length to PLS1
pls1.genes = subset(pls1.genes, !is.na(pls1.genes$z_uncorr))
pls1.genes = pls1.genes[order(pls1.genes$z_uncorr, decreasing=T),]

genes= pls1.genes[,c(1:7)]    # Background gene list (i.e. includes gene length)
pls1 = pls1.genes[,c(1,7,8)]

# # keep genes in corr only
lake2 = lake
lake2$hgnc_symbol = lake$Gene
lake2 = lake2 %>% left_join(pls1, by='hgnc_symbol')
lake2 = subset(lake2, !is.na(lake2$z_uncorr))

# create null models
script_dir = "C:\\Users\\lihon\\Desktop\\abcd_study\\pnc\\gene_decoding\\"
source(paste0(script_dir, "\\gene_nulls_nearest_neighbour.R")) # Gene enrichment function
source(paste0(script_dir, "\\run_gene_enrichment.R"))          # Wrapper script to submit gene lists to gene enrichment function
source(paste0(script_dir, "\\plot_enrichment_ggridges.R"))     # Plot enrichment function

n_perm = 5000 # Number of random permutations
out.path = "C:\\Users\\lihon\\Desktop\\abcd_study\\pnc\\vis_v2_r400_ro1\\corr_gene_enrichment_cell\\"

################ Adult (Post Natal) Cell Type Enrichment ################
postnatal_res = run_gene_enrichment(lake2, pls1, TRUE, genes, n_perm, out.path)
p.postnatal<-plot_enrichment_ggridges(postnatal_res$gene_stats, postnatal_res$gene_null, TRUE, -1000, 1000) 
pdf(paste0(out.path, "postnatal.enrichment.pdf"), 4, 4); p.postnatal; dev.off()
print(xtable(postnatal_res$gene_length_stats, booktabs=TRUE), include.rownames=FALSE)


