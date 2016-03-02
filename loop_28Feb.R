rm(list=ls())
#setwd("C:/Katies Data/Amphibians/P1_phylog_synthesis")
Marske_phylo <- function(species_vector, raw_data_path, output_path){
  library(ape)  
  library(pegas)
  library(adegenet)
  library(strataG)
  library(sidier)
  library(reshape2)
  # Creates an empty data frame to save the species genetic stats
  Marske_stats <- as.data.frame(matrix(nrow=length(species_vector), ncol=10))
  colnames(Marske_stats) <- c("Species", "N_seqs", "N_haps", "nuc_div", "Var_nuc_div", "TajD", "TajD.p.beta", "R2", "R2.p", "Phist")
  for (sp in seq_along(species_vector)){
    #gets seq data in as DNAbin
    setwd(paste(raw_data_path, "/Fasta/", sep=""))
    sp_align_temp <-read.dna(paste(species_vector[sp],"_align.fasta", sep=""), format = "fasta", as.character = FALSE, as.matrix=NULL)
    #reads csv file for each species
    setwd(paste(raw_data_path, "/Csv/", sep=""))
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
    sp_nuc.div_temp  <-  nuc.div(sp_align_temp, variance = T, pairwise.deletion = FALSE)
    #sp_hap.div<-haplotypic.diversity(sp_gtype) #haplotype diversity, strataG
    #Tajima's D, pegas, #Pval.beta - p-value assuming that D follows a beta distribution (Tajima, 1989)
    sp_tajD_temp <-tajima.test(sp_align_temp) 
    # R squared test, Ramos-Onsins-Rozas Test of Neutrality (Pegas)
    sp_R2_temp   <-R2.test(sp_align_temp, B=10000, plot=FALSE, quiet=TRUE)
    
    ### Distance matrix ####
    sp_gendist_temp  <-  dist.dna(sp_align_temp, model = "K80", variance = FALSE, gamma = FALSE, pairwise.deletion = FALSE, base.freq = NULL, as.matrix = FALSE)
    #using reshape2
    pairwise.gendist_temp <-melt(as.matrix(sp_gendist_temp), varnames = c("row", "col"), value.name="gendist")
    setwd(paste(output_path, "/Gendist/", sep=""))
    write.csv(pairwise.gendist_temp, paste(species_vector[sp],"_gendist.csv", sep=""), row.names=F)
    
    ### strataG ####
    #save haps in a DNAbin
    setwd(paste(output_path, "/Hap_fasta/", sep=""))
    sp_haps_temp <- GetHaplo(align=sp_align_temp, saveFile=T, outname=paste(species_vector[sp], "_haps.fasta", sep=""), format="fasta", seqsNames="Inf.Hap") 
    sp_ghaps_temp <-read.fasta(paste(species_vector[sp], "_haps.fasta", sep=""))
    #contents of row and col columns are UniqueInds
    #need the sep element to not put space in the name
    Marske_stats[sp,] <- c(species_vector[sp], N_seqs_temp, N_haps_temp, sp_nuc.div_temp[1], sp_nuc.div_temp[2], sp_tajD_temp$D, sp_tajD_temp$Pval.beta, sp_R2_temp$R2, sp_R2_temp$P.val, NA)
  }
setwd(output_path)
write.csv(Marske_stats, "results.csv")
}

### A function to estimate all the parameters you need ####
## the parameters are: 1) a list of species, 2) the path to your raw (fasta and csv) data and 3) path to the folder to store the results
# to run the function
# makes the species vector, you can put any vector of species (e.g. just one: species_vector <- "Ambystoma_bishopi")
setwd("/Users/afr/Desktop/Marske_2016/Marske_clean//Raw_data/Csv/")
# to include all the species in the folder of raw data csv
species_vector <- list.files(pattern = ".csv", full.names=FALSE)
species_vector <- gsub(".csv", "", species_vector)
# defines the path for the raw data, this path is for the folder containing the folders "fasta" and "csv"
raw_data_path <- "~/Desktop/Marske_2016/Marske_clean/Raw_data"
# defines the path for the output, this path is for the folder containing the folders "gendist" and "hap_fasta"
output_path <- "~/Desktop/Marske_2016/Marske_clean/R_outputs/"

Marske_phylo(species_vector, raw_data_path, output_path)



