rm(list=ls())
#setwd("C:/Katies Data/Amphibians/P1_phylog_synthesis")
<<<<<<< HEAD
Marske_phylo <- function(species_vector, raw_data_path, output_path){
=======
Marske_phylo <- function(species_vector){
>>>>>>> 22bf18a71fc34f2b56cab935b7203bfffaeb1158
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
<<<<<<< HEAD
    setwd(paste(raw_data_path, "/Fasta/", sep=""))
    sp_align_temp <-read.dna(paste(species_vector[sp],"_align.fasta", sep=""), format = "fasta", as.character = FALSE, as.matrix=NULL)
    #reads csv file for each species
    setwd(paste(raw_data_path, "/Csv/", sep=""))
=======
    setwd("/Users/afr/Desktop/Marske_2016/P1_phylog_synthesis/Raw_data/Fasta/")
    sp_align_temp <-read.dna(paste(species_vector[sp],"_align.fasta", sep=""), format = "fasta", as.character = FALSE, as.matrix=NULL)
    #reads csv file for each species
    setwd("/Users/afr/Desktop/Marske_2016/P1_phylog_synthesis/Raw_data/Csv/")
>>>>>>> 22bf18a71fc34f2b56cab935b7203bfffaeb1158
    sp_data_temp  <-read.csv(paste(species_vector[sp],".csv", sep=""), stringsAsFactors = F)
    #Sorts sp_data in the order of sp_align
    sp_data_temp <- sp_data_temp[match(rownames(sp_align_temp), sp_data_temp$UniqueInd),] 
    
    ### Haplotype stuff ####
<<<<<<< HEAD
    #get unique haps using sidier and save it in the outputs directory  
=======
    #get unique haps using sidier and save it in the outputs directory
>>>>>>> 22bf18a71fc34f2b56cab935b7203bfffaeb1158
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
<<<<<<< HEAD
    setwd(paste(output_path, "/Gendist/", sep=""))
=======
    setwd("/Users/afr/Desktop/Marske_2016/P1_phylog_synthesis/R_outputs/Gendist/")
>>>>>>> 22bf18a71fc34f2b56cab935b7203bfffaeb1158
    write.csv(pairwise.gendist_temp, paste(species_vector[sp],"_gendist.csv", sep=""), row.names=F)
    
    ### strataG ####
    #save haps in a DNAbin
<<<<<<< HEAD
    setwd(paste(output_path, "/Hap_fasta/", sep=""))
=======
    setwd("/Users/afr/Desktop/Marske_2016/P1_phylog_synthesis/R_outputs/Hap_fasta/")
>>>>>>> 22bf18a71fc34f2b56cab935b7203bfffaeb1158
    sp_haps_temp <- GetHaplo(align=sp_align_temp, saveFile=T, outname=paste(species_vector[sp], "_haps.fasta", sep=""), format="fasta", seqsNames="Inf.Hap") 
    sp_ghaps_temp <-read.fasta(paste(species_vector[sp], "_haps.fasta", sep=""))
    #contents of row and col columns are UniqueInds
    #need the sep element to not put space in the name
    Marske_stats[sp,] <- c(species_vector[sp], N_seqs_temp, N_haps_temp, sp_nuc.div_temp, sp_tajD_temp$D, sp_tajD_temp$Pval.beta, sp_R2_temp$R2, sp_R2_temp$P.val, "phist")
  }
<<<<<<< HEAD
setwd(output_path)
=======
>>>>>>> 22bf18a71fc34f2b56cab935b7203bfffaeb1158
write.csv(Marske_stats, "results.csv")
}

### A function to estimate all the parameters you need ####
<<<<<<< HEAD
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
=======
## the parameters are: 1) a list of species, 2) a list of parameters to estimate (e.g. phist, TD, R2...)
# makes the species vector
setwd("/Users/afr/Desktop/Marske_2016/P1_phylog_synthesis/Raw_data/Csv/")
species_vector <- list.files(pattern = ".csv", full.names=FALSE)
species_vector <- gsub(".csv", "", species_vector)
setwd("/Users/afr/Desktop/Marske_2016/P1_phylog_synthesis/R_outputs/")
Marske_phylo(species_vector[-2])
>>>>>>> 22bf18a71fc34f2b56cab935b7203bfffaeb1158



