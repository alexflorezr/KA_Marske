rm(list=ls())
#setwd("C:/Katies Data/Amphibians/P1_phylog_synthesis")
Marske_phylo <- function(species_vector){
  library(ape)  
  library(pegas)
  library(adegenet)
  library(strataG)
  library(sidier)
  library(reshape2)
  # Creates an empty data frame to save the species genetic stats
  Marske_stats <- as.data.frame(matrix(nrow=length(species_vector), ncol=9))
  colnames(Marske_stats) <- c("Species", "N_seqs", "N_haps", "nuc_div", "TajD", "TajD.p.beta", "R2", "R2.p", "Phist")
  for (sp in seq_along(species_vector)){
    #gets seq data in as DNAbin
    setwd("/Users/afr/Desktop/Marske_2016/P1_phylog_synthesis/Raw_data/Fasta/")
    sp_align_temp <-read.dna(paste(species_vector[sp],"_align.fasta", sep=""), format = "fasta", as.character = FALSE, as.matrix=NULL)
    #reads csv file for each species
    setwd("/Users/afr/Desktop/Marske_2016/P1_phylog_synthesis/Raw_data/Csv/")
    sp_data_temp  <-read.csv(paste(species_vector[sp],".csv", sep=""), stringsAsFactors = F)
    #Sorts sp_data in the order of sp_align
    sp_data_temp <- sp_data_temp[match(rownames(sp_align_temp), sp_data_temp$UniqueInd),] 
    
    ### Haplotype stuff ####
    #get unique haps using sidier and save it in the outputs directory
    haps_temp <-FindHaplo(align=sp_align_temp, saveFile=F)
    colnames(haps_temp)<-c("UniqueInd", "haplotype")
    full_data_temp <-merge(sp_data_temp, haps_temp)
    
    ### Diversity stats ####
    # Number of sequences
    N_seqs_temp  <- length(sp_data_temp$UniqueInd)
    # Number of haplotypes
    N_haps_temp  <- length(unique(full_data_temp$haplotype))
    # Nucleotide diversity (Pegas)
    sp_nuc.div_temp  <-  nuc.div(sp_align_temp, variance = F, pairwise.deletion = FALSE)
    #sp_hap.div<-haplotypic.diversity(sp_gtype) #haplotype diversity, strataG
    #Tajima's D, pegas, #Pval.beta - p-value assuming that D follows a beta distribution (Tajima, 1989)
    sp_tajD_temp <-tajima.test(sp_align_temp) 
    # R squared test, Ramos-Onsins-Rozas Test of Neutrality (Pegas)
    sp_R2_temp   <-R2.test(sp_align_temp, B=10000, plot=FALSE, quiet=TRUE)
    
    ### Distance matrix ####
    sp_gendist_temp  <-  dist.dna(sp_align_temp, model = "K80", variance = FALSE, gamma = FALSE, pairwise.deletion = FALSE, base.freq = NULL, as.matrix = FALSE)
    #using reshape2
    pairwise.gendist_temp <-melt(as.matrix(sp_gendist_temp), varnames = c("row", "col"), value.name="gendist")
    setwd("/Users/afr/Desktop/Marske_2016/P1_phylog_synthesis/R_outputs/Gendist/")
    write.csv(pairwise.gendist_temp, paste(species_vector[sp],"_gendist.csv", sep=""), row.names=F)
    
    ### strataG ####
    #save haps in a DNAbin
    setwd("/Users/afr/Desktop/Marske_2016/P1_phylog_synthesis/R_outputs/Hap_fasta/")
    sp_haps_temp <- GetHaplo(align=sp_align_temp, saveFile=T, outname=paste(species_vector[sp], "_haps.fasta", sep=""), format="fasta", seqsNames="Inf.Hap") 
    sp_ghaps_temp <-read.fasta(paste(species_vector[sp], "_haps.fasta", sep=""))
    #contents of row and col columns are UniqueInds
    #need the sep element to not put space in the name
    Marske_stats[sp,] <- c(species_vector[sp], N_seqs_temp, N_haps_temp, sp_nuc.div_temp, sp_tajD_temp$D, sp_tajD_temp$Pval.beta, sp_R2_temp$R2, sp_R2_temp$P.val, "phist")
  }
write.csv(Marske_stats, "results.csv")
}

### A function to estimate all the parameters you need ####
## the parameters are: 1) a list of species, 2) a list of parameters to estimate (e.g. phist, TD, R2...)
# makes the species vector
setwd("/Users/afr/Desktop/Marske_2016/P1_phylog_synthesis/Raw_data/Csv/")
species_vector <- list.files(pattern = ".csv", full.names=FALSE)
species_vector <- gsub(".csv", "", species_vector)
setwd("/Users/afr/Desktop/Marske_2016/P1_phylog_synthesis/R_outputs/")
Marske_phylo(species_vector[-2])



