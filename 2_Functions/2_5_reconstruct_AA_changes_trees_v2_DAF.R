########################################################################################################################################
## Reconstruct AA changes in trees
########################################################################################################################################
library(stringr)
library(ape)
library(phytools)
library(coda)
library(vcfR)
library(lubridate)
library(doParallel)
library(phangorn)
library(seqinr)
library(here)
library(roxygen2)
library(docstring)
library(tidyverse)
library(readxl)



########################################################################################################################################
## Load data
########################################################################################################################################
here::i_am("DENV_3_fitness.Rmd")

load("6_Session_Data/Rebuild03_18_25.Rdata")


### DENV 1 ####
## Input fasta data (PHYDAT format)
tree_data_DENV1 <- read.phyDat("1_Data/1_5_DENV/Seq_vcf_data/denv_seqs_wgs/whole_genome/d1_n1026.fas",
                              format = "fasta")
names(tree_data_DENV1) <- gsub(" ", "_", names(tree_data_DENV1))

## Input fasta data (SEQINR format)
fasta_data_DENV1 = read.fasta('1_Data/1_5_DENV/Seq_vcf_data/denv_seqs_wgs/whole_genome/d1_n1026.fas', forceDNAtolower = T, whole.header = T)
data_seq_DENV1 = getSequence(fasta_data_DENV1)
data_name_seq_DENV1 = names(fasta_data_DENV1)
data_name_seq_DENV1 = gsub(" ", "_", data_name_seq_DENV1)


tree_1 = tree_DENV1

### DENV 2 ####
tree_data_DENV2 <- read.phyDat("1_Data/1_5_DENV/Seq_vcf_data/denv_seqs_wgs/whole_genome/d2_n796.fas",
                               format = "fasta")
names(tree_data_DENV2) <- gsub(" ", "_", names(tree_data_DENV2))

## Input fasta data (SEQINR format)
fasta_data_DENV2 = read.fasta('1_Data/1_5_DENV/Seq_vcf_data/denv_seqs_wgs/whole_genome/d2_n796.fas', forceDNAtolower = T, whole.header = T)
data_seq_DENV2 = getSequence(fasta_data_DENV2)
data_name_seq_DENV2 = names(fasta_data_DENV2)
data_name_seq_DENV2 = gsub(" ", "_", data_name_seq_DENV2)


tree_2 = tree_DENV2

### DENV 3 ####
tree_data_DENV3 <- read.phyDat("1_Data/1_5_DENV/Seq_vcf_data/denv_seqs_wgs/whole_genome/d3_n625_degap.fas",
                               format = "fasta")
names(tree_data_DENV3) <- gsub(" ", "_", names(tree_data_DENV3))

## Input fasta data (SEQINR format)
fasta_data_DENV3 = read.fasta("1_Data/1_5_DENV/Seq_vcf_data/denv_seqs_wgs/whole_genome/d3_n625_degap.fas", forceDNAtolower = T, whole.header = T)
data_seq_DENV3 = getSequence(fasta_data_DENV3)
data_name_seq_DENV3 = names(fasta_data_DENV3)
data_name_seq_DENV3 = gsub(" ", "_", data_name_seq_DENV3)


tree_3 = tree_DENV3


### DENV 4 ####
tree_data_DENV4 <- read.phyDat("1_Data/1_5_DENV/Seq_vcf_data/denv_seqs_wgs/whole_genome/d4_n477.fas",
                               format = "fasta")
names(tree_data_DENV4) <- gsub(" ", "_", names(tree_data_DENV4))

## Input fasta data (SEQINR format)
fasta_data_DENV4 = read.fasta("1_Data/1_5_DENV/Seq_vcf_data/denv_seqs_wgs/whole_genome/d4_n477.fas", forceDNAtolower = T, whole.header = T)
data_seq_DENV4 = getSequence(fasta_data_DENV4)
data_name_seq_DENV4 = names(fasta_data_DENV4)
data_name_seq_DENV4 = gsub(" ", "_", data_name_seq_DENV4)


tree_4 = tree_DENV4

# OR IF trees aren't already loaded in ####
# #Read in tree as phylo 
# tree_DENV1 = read.nexus("1_Data/1_5_DENV/best_tree_d1_Lin.nexus")
# tree_DENV1 = collapse.singles(ladderize(multi2di(tree_phylo_DENV1, random = F), right = F))

#setdiff( names(tree_data_DENV1), tree$tip.label)
# 934 tip labels vs 1026 in the treedata


########################################################################################################################################

########################################################################################################################################
## Prepare data: make the IDs match
########################################################################################################################################]


### DENV 1 ####
correspondence_ID_EPI_genious_DENV1 = read.csv("1_Data/1_5_DENV/Seq_vcf_data/denv_seqs_wgs/whole_genome/d1_n1026_geneious_metadata.csv", 
                                               col.names = c("Name",	"Created Date",	"Description",	"Imported From: Filename",	"Imported From: Path",	"Sequence", "A", "B"))
 

# Remove the dates with multiples designation in metadata names  ("2012.361 2" → "2012.361") )
correspondence_ID_EPI_genious_DENV1$Name <- gsub("^([0-9]+\\.[0-9]+).*", "\\1", correspondence_ID_EPI_genious_DENV1$Name) 
  # ^([0-9]+\\.[0-9]+).*
  #  ^ → Start of string.
  #  ([0-9]+\\.[0-9]+) → Captures a pattern of a year followed by a decimal (e.g., "2011.501").
  # .* → Matches anything that follows.
  #  "\\1" → Keeps only the first matched group (the year and decimal part), removing everything else.



### DENV 2 ####
correspondence_ID_EPI_genious_DENV2 = read.csv("1_Data/1_5_DENV/Seq_vcf_data/denv_seqs_wgs/whole_genome/d2_n796_geneious_metadata.csv", 
                                               col.names = c("Name",	"Created Date",	"Description",	"Imported From: Filename",	"Imported From: Path",	"Sequence", "A", "B"))


# Remove the dates with multiples designation in metadata names  ("2012.361 2" → "2012.361") )
correspondence_ID_EPI_genious_DENV2$Name <- gsub("^([0-9]+\\.[0-9]+).*", "\\1", correspondence_ID_EPI_genious_DENV2$Name) 
# ^([0-9]+\\.[0-9]+).*
#  ^ → Start of string.
#  ([0-9]+\\.[0-9]+) → Captures a pattern of a year followed by a decimal (e.g., "2011.501").
# .* → Matches anything that follows.
#  "\\1" → Keeps only the first matched group (the year and decimal part), removing everything else.

### DENV 3 ####
correspondence_ID_EPI_genious_DENV3 = read.csv("1_Data/1_5_DENV/Seq_vcf_data/denv_seqs_wgs/whole_genome/d3_n625_degap_geneious_metadata.csv", 
                                               col.names = c("Name",	"Created Date",	"Description",	"Imported From: Filename",	"Imported From: Path",	"Sequence", "A", "B"))


# Remove the dates with multiples designation in metadata names  ("2012.361 2" → "2012.361") )
correspondence_ID_EPI_genious_DENV3$Name <- gsub("^([0-9]+\\.[0-9]+).*", "\\1", correspondence_ID_EPI_genious_DENV3$Name) 
# ^([0-9]+\\.[0-9]+).*
#  ^ → Start of string.
#  ([0-9]+\\.[0-9]+) → Captures a pattern of a year followed by a decimal (e.g., "2011.501").
# .* → Matches anything that follows.
#  "\\1" → Keeps only the first matched group (the year and decimal part), removing everything else.


### DENV 4 ####
correspondence_ID_EPI_genious_DENV4 = read.csv("1_Data/1_5_DENV/Seq_vcf_data/denv_seqs_wgs/whole_genome/d4_n477_geneious_metadata.csv", 
                                               col.names = c("Name",	"Created Date",	"Description",	"Imported From: Filename",	"Imported From: Path",	"Sequence", "A", "B"))


# Remove the dates with multiples designation in metadata names  ("2012.361 2" → "2012.361") )
correspondence_ID_EPI_genious_DENV4$Name <- gsub("^([0-9]+\\.[0-9]+).*", "\\1", correspondence_ID_EPI_genious_DENV4$Name) 
# ^([0-9]+\\.[0-9]+).*
#  ^ → Start of string.
#  ([0-9]+\\.[0-9]+) → Captures a pattern of a year followed by a decimal (e.g., "2011.501").
# .* → Matches anything that follows.
#  "\\1" → Keeps only the first matched group (the year and decimal part), removing everything else.


correspondence_ID_EPI_all <- list(correspondence_ID_EPI_genious_DENV1, 
                                  correspondence_ID_EPI_genious_DENV2, 
                                  correspondence_ID_EPI_genious_DENV3, 
                                  correspondence_ID_EPI_genious_DENV4 )

names(correspondence_ID_EPI_all) <- c("DENV1", "DENV2", "DENV3", "DENV4")

# Apply transformations to all data frames in the list, create longname and make it look liek the vcf file names
correspondence_ID_EPI_all <- lapply(correspondence_ID_EPI_all, function(df) {
  # Create the 'longname' column dynamically
  df$longname <- gsub(" ", "_", paste(df$Name, df$Description, sep = "_"))
  
  # Reorder columns to put 'longname' first
  df <- df[, c("longname", setdiff(names(df), "longname"))]
  
  return(df)
})



########################################################################################################################################

########################################################################################################################################
## DENV Reconstruct ancestral sequences within tree
########################################################################################################################################
process_pml_trees <- function(i_range = 1:4) {
  #' Process Multiple Trees with Maximum Likelihood and Ancestral Reconstruction
  #'
  #' @description This function iterates over a range of tree indices, fitting a likelihood tree model (GTR) and performing Bayesian ancestral reconstruction.
  #'
  #' @param i_range Numeric vector. A sequence of indices representing tree versions to process (default: `1:4`).
  #'
  #' @return A list containing:
  #' \itemize{
  #'   \item{\code{fits}}{A list of fitted PML objects for each tree.}
  #'   \item{\code{ancestral_states}}{A list of Bayesian ancestral reconstructions for each tree.}
  #' }
  #'
  #' @examples
  #' result <- process_pml_trees(1:4)
  #' fit_1 <- result$fits[[1]]
  #' anc.bayes_1 <- result$ancestral_states[[1]]
  #' 
  #' @export
  
  library(phangorn)
  
  # Pre-allocate lists
  fits <- vector("list", length(i_range))
  ancestral_states <- vector("list", length(i_range))
  
  # Iterate over trees
  for (i in i_range) {
    tree_var <- paste0("tree_", i)
    data_var <- paste0("tree_data_DENV", i)
    
    if (!exists(tree_var) || !exists(data_var)) {
      print(paste("Skipping i =", i, "- Missing tree or data variable"))
      next
    }
    
    tree <- get(tree_var)
    data <- get(data_var)
    
    print(paste("Processing tree", i, "with data:", data_var))
    
    # Fit the likelihood tree model (GTR)
    fit <- pml(tree = tree, data = data)
    fit <- optim.pml(fit, model = "GTR", control = pml.control(trace = 0))
    
    # Bayesian ancestral state reconstruction
    anc_bayes <- ancestral.pml(fit, "bayes")
    
    # Store results
    fits[[i]] <- fit
    ancestral_states[[i]] <- anc_bayes
  }
  
  # Return structured output
  return(list(fits = fits, ancestral_states = ancestral_states))
}

pml_results <- process_pml_trees(1:4)
names(pml_results$fits) <- c("DENV1", "DENV2", "DENV3", "DENV4")
names(pml_results$ancestral_states) <- c("DENV1", "DENV2", "DENV3", "DENV4")

fit_1 <- pml_results$fits[[1]]
anc.bayes_1 <- pml_results$ancestral_states[[1]]

fit_2 <- pml_results$fits[[2]]
anc.bayes_2 <- pml_results$ancestral_states[[2]]

fit_3 <- pml_results$fits[[3]]
anc.bayes_3 <- pml_results$ancestral_states[[3]]

fit_4 <- pml_results$fits[[4]]
anc.bayes_4 <- pml_results$ancestral_states[[4]]


# Old Code for reconstruction  --------------------------------------------

# tree_1
# 
# # pml() A phangorn package that create a likelyhood of a tree, requires a phylo class tree and then a phydat data object
# # optim.pml allows using a GTR model for evolution
# fit_1 <- pml(tree = tree_1, data = tree_data_DENV1)
# fit_1 <- optim.pml(fit_1, model="GTR", control = pml.control(trace=0))
# 
# setdiff( names(tree_data_DENV1), names(fit_1$data))
# 
# 
# 
# # Marginal reconstruction of the ancestral character states.
# anc.bayes_1 <- ancestral.pml(fit_1, "bayes")






########################################################################################################################################

########################################################################################################################################
## DENV Write down matrix of sequences (including nodes)
########################################################################################################################################
process_matrix_data_bin <- function(i_range = 1:4) {
  #' Process Bayesian Ancestral State Reconstruction Matrices
  #'
  #' @description This function iterates over multiple Bayesian ancestral reconstruction datasets (`anc.bayes_i`)
  #' and creates corresponding `matrix_data_bin_i` and `matrix_data_i` matrices.
  #'
  #' @param i_range Numeric vector. A sequence of indices representing datasets to process (default: `1:4`).
  #'
  #' @return A list containing:
  #' \itemize{
  #'   \item{\code{matrix_data_bin}}{A list of binary state matrices for each dataset.}
  #'   \item{\code{matrix_data}}{A list of character matrices for reconstructed nucleotides.}
  #' }
  #'
  #' @examples
  #' result <- process_matrix_data_bin(1:4)
  #' matrix_data_bin_1 <- result$matrix_data_bin[[1]]
  #' matrix_data_1 <- result$matrix_data[[1]]
  #'
  #' @export
  
  library(phangorn)
  
  # Pre-allocate lists to store results
  matrix_data_bin_list <- vector("list", length(i_range))
  matrix_data_list <- vector("list", length(i_range))
  
  # Iterate over datasets
  for (i in i_range) {
    anc_var <- paste0("anc.bayes_", i)
    seq_var <- paste0("data_seq_DENV", i)
    
    if (!exists(anc_var) || !exists(seq_var)) {
      print(paste("Skipping i =", i, "- Missing ancestral or sequence data"))
      next
    }
    
    anc_bayes <- get(anc_var)
    data_seq <- get(seq_var)
    
    print(paste("Processing dataset", i, "with ancestral data:", anc_var))
    
    # Initialize matrices
    matrix_data_bin <- matrix_data <- matrix(0, length(names(anc_bayes)), length(data_seq[[1]]))
    
    # Process each sequence
    for (j in seq_len(length(names(anc_bayes)))) {
      matrix_data_bin[j, ] <- apply(anc_bayes[[j]], MAR = 1, function(x) {
        if (all(is.na(x)) || sum(x, na.rm = TRUE) == 0) {
          return(NA)  # Avoid errors by returning NA when necessary
        }
        x <- x / sum(x, na.rm = TRUE)  # Normalize probabilities
        tmp <- max(x, na.rm = TRUE)  # Get max probability
        tmp_which <- which.max(x)  # Index of most probable nucleotide
        
        if (!is.na(tmp) && tmp > 0.9) {
          return(tmp_which)
        } else {
          return(NA)
        }
      })[attr(anc_bayes, "index")]
      
      tmp <- factor(matrix_data_bin[j,], levels = 1:4)
      levels(tmp) <- attr(anc_bayes, 'levels')
      matrix_data[j, ] <- as.character(tmp)
      matrix_data[j, which(is.na(matrix_data[j, ]))] <- 'n'  # Assign 'n' to missing values
    }
    
    # Assign row names
    rownames(matrix_data) <- names(anc_bayes)
    
    # Store results
    matrix_data_bin_list[[i]] <- matrix_data_bin
    matrix_data_list[[i]] <- matrix_data
  }
  
  # Return structured output
  return(list(matrix_data_bin = matrix_data_bin_list, matrix_data = matrix_data_list))
}


matrix_data_results <- process_matrix_data_bin(1:4)
names(matrix_data_results$matrix_data) <- c("DENV1", "DENV2", "DENV3", "DENV4")
names(matrix_data_results$matrix_data_bin) <- c("DENV1", "DENV2", "DENV3", "DENV4")

# Access results for dataset 1
matrix_data_bin_1 <- matrix_data_results$matrix_data_bin$DENV1
matrix_data_1 <- matrix_data_results$matrix_data$DENV1

# Access results for dataset 2, 3, etc
matrix_data_bin_2 <- matrix_data_results$matrix_data_bin$DENV2
matrix_data_2 <- matrix_data_results$matrix_data$DENV2

matrix_data_bin_3 <- matrix_data_results$matrix_data_bin$DENV3
matrix_data_3 <- matrix_data_results$matrix_data$DENV3

matrix_data_bin_4 <- matrix_data_results$matrix_data_bin$DENV4
matrix_data_4 <- matrix_data_results$matrix_data$DENV4

# Old Code for Matrix -----------------------------------------------------
# 
# # Creates two matrices, matrix_data_bin and matrix_data, both initialized with zeros.
# matrix_data_bin_1 = matrix_data_1 = matrix(0, length(names(anc.bayes_1)), length(data_seq_DENV1[[1]]))
# 
# for(i in 1:length(names(anc.bayes_1))) {
#   matrix_data_bin_1[i,] = apply(anc.bayes_1[[i]], MAR = 1, function(x) {
#     if (all(is.na(x)) || sum(x, na.rm = TRUE) == 0) {
#       return(NA)  # If all values are NA or sum is zero, return NA to avoid errors
#     }
#     x = x / sum(x, na.rm = TRUE)  # Normalize probabilities safely
#     tmp = max(x, na.rm = TRUE)  # Find the max probability safely
#     tmp_which = which.max(x)  # Get the index of the most probable nucleotide
#     
#     if (!is.na(tmp) && tmp > 0.9) {
#       return(tmp_which)
#     } else {
#       return(NA)
#     }
#   })[attr(anc.bayes_1, "index")]
#   
#   tmp = factor(matrix_data_bin_1[i,], levels = 1:4)
#   levels(tmp) = attr(anc.bayes_1, 'levels')
#   matrix_data_1[i,] = as.character(tmp)
#   matrix_data_1[i, which(is.na(matrix_data_1[i,]))] = 'n' #Assigns n if the value is NA
# }
# 
# 
# 
# rownames(matrix_data_1) = names(anc.bayes_1)
# View(head(matrix_data_1,100))

# remove(tree_DENV1) 
# remove(fasta_data_DENV1) 
# remove(data_seq_DENV1)
########################################################################################################################################

########################################################################################################################################
## Check for 1 position that everything is fine, this can take a minute to render
# NS2A is highly variable so check there? 3480:4133, set a SNP position column to look at in matrix 
########################################################################################################################################
# Get the number of unique values per column (this takes some time)

for( i in 1:4) {
  matrix_data <- get(paste0("matrix_data_bin_", i))
  
  position_variability_df <- data.frame(
    position = seq_len(ncol(matrix_data)),       # 1 to number of columns
    unique_bases = apply(matrix_data, 2, function(col) length(unique(col)))
  )
  
  high_variance_pos <- which(position_variability_df$unique_bases %in% c(4, 5))
  
  max_unique_col <- which.max(position_variability_df$unique_bases)
  
  # --- Histogram Plot ---
  png_filename <- paste0("DENV", i, "_unique_base_histogram.png")
  png(here("3_Output_Figures", png_filename), width = 5, height = 5, units = "in", res = 300)
  
  hist_obj <- hist(position_variability_df$unique_bases,
                   main = paste("DENV", i, "Histogram of Unique Base Counts"),
                   xlab = "Number of Unique Bases",
                   ylab = "Frequency",
                   col = "lightblue",
                   border = "black")
  
  nonzero <- hist_obj$counts > 0
  text(x = hist_obj$mids[nonzero],
       y = hist_obj$counts[nonzero],
       labels = hist_obj$counts[nonzero],
       pos = 3, cex = 0.8, col = "black")
  
  dev.off()  # Save PNG
  
  
  
  p <- plot_ly(data = position_variability_df,
          x = ~position,
          y = ~unique_bases,
          type = 'scatter',
          mode = 'markers',
          marker = list(color = 'darkgreen', size = 3),
          text = ~paste("Position:", position, "<br>Unique Bases:", unique_bases),
          hoverinfo = 'text') %>%
    layout(title = paste0("DENV",i," Histogram of Unique Base Counts per Column"),
           xaxis = list(title = "Genome Position (Column Index)"),
           yaxis = list(title = "Number of Unique Bases"))
  
  htmlwidgets::saveWidget(as_widget(p), here("3_Output_Figures", paste0("DENV",i," Histogram of Unique Base Counts per Column", ".html")))
}


#This plot takes forever
# plot.phylo(tree_DENV1, show.tip.label = F, node.color = matrix_data_bin_1[,231])


########################################################################################################################################

########################################################################################################################################
## Split this matrix into ORFs, split by codons and translate to AA (### DAF Optimized code)
########################################################################################################################################


# Old code ----------------------------------------------------------------


# # Load ORF position data from a CSV file
# data_orfs = read.csv("1_Data/1_5_DENV/Position_genes_Dengue.csv")
# 
# # Pre-allocate lists
# list_matrices_ORFs_AA = vector("list", length = nrow(data_orfs))
# list_matrices_ORFs_codon = vector("list", length = nrow(data_orfs))
# 
# # Iterate over ORFs
# for (i in seq_len(nrow(data_orfs))) {
#   print(paste("Processing ORF", i, "of", nrow(data_orfs)))  # Track progress
#   start_idx <- max(1, data_orfs$start[i])  
#   end_idx <- min(ncol(matrix_data_1), data_orfs$end[i])
#   
#   # Extract ORF data once to avoid repeated indexing
#   orf_data <- matrix_data_1[, start_idx:end_idx]
#   
#   # AA sequences (optimized apply usage)
#   list_matrices_ORFs_AA_1[[i]] <- t(apply(orf_data, 1, translate))
#   
#   # Precompute codon matrix size
#   num_codons <- (end_idx - start_idx + 1) / 3
#   list_matrices_ORFs_codon_1[[i]] <- matrix(NA, nrow = nrow(matrix_data_1), ncol = num_codons)
#   
#   # Efficient codon extraction using matrix slicing
#   for (j in seq_len(num_codons)) {
#     col_indices <- (1:3) + (j - 1) * 3
#     list_matrices_ORFs_codon_1[[i]][, j] <- apply(orf_data[, col_indices], 1, paste0, collapse = "")
#   }
#   
#   # Assign row names once (avoid redundant operations)
#   rownames(list_matrices_ORFs_AA_1[[i]]) <- names(anc.bayes_1)
#   rownames(list_matrices_ORFs_codon_1[[i]]) <- names(anc.bayes_1)
# }
########


process_ORFs_into_AA_and_Codons <- function(data_orfs, matrix_data, anc_bayes, orf_index = NULL, offset = 0) {
  #' Process a matrix of sequence data into codons and amino acids
  #'
  #' @description This function processes a matrix of nucleotide sequence data into codon and amino acid matrices.
  #' It groups sequences based on open reading frames (ORFs) or genes, as specified in a genome map provided by the user.
  #'
  #' @param data_orfs Character or data.frame. A path to a CSV file containing ORF positions or a preloaded data.frame.
  #' The ORF file should contain columns `start` and `end` that define the ORF locations.
  #' @param matrix_data Matrix. A matrix of nucleotide sequences, where rows represent different sequences and columns represent nucleotide positions.
  #' @param anc_bayes Named vector. A reference sequence or ancestral state vector used to assign row names to output matrices.
  #' @param orf_index Integer or NULL. If provided, only the specified ORF index will be processed. If NULL (default), all ORFs are processed.
  #' @param offset Moves the start and end of genes/orfs to align with sequences.
  #' 
  #' @return A list containing:
  #' \itemize{
  #'   \item{`AA_matrices`}{A list of matrices, where each matrix contains translated amino acid sequences for an ORF.}
  #'   \item{`codon_matrices`}{A list of matrices, where each matrix contains codon sequences extracted from `matrix_data`.}
  #' }
  #'
  #' @details ORF regions are truncated to the nearest multiple of three to ensure correct codon alignment.
  #'
  #' @examples
  #' output <- process_ORFs_into_AA_and_Codons(data_orfs, matrix_data, anc_bayes)
  #' list_matrices_ORFs_AA_1 <- output$AA_matrices
  #' list_matrices_ORFs_codon_1 <- output$codon_matrices
  #'
  
  # Read ORF position data if file path is provided
  if (is.character(data_orfs)) {
    data_orfs <- read.csv(data_orfs)
  }
  
  # Determine the ORF indices to process
  if (is.null(orf_index)) {
    orf_indices <- seq_len(nrow(data_orfs))
  } else {
    orf_indices <- orf_index
  }
  
  # Pre-allocate lists
  AA_matrices <- vector("list", length = nrow(data_orfs))
  codon_matrices <- vector("list", length = nrow(data_orfs))
  
  # Iterate over ORFs
  for (i in orf_indices) {
    print(paste("Processing ORF", i, "of", nrow(data_orfs)))
    
    start_idx <- max(1, data_orfs$Genome_Start[i] + offset)
    end_idx <- min(ncol(matrix_data), data_orfs$Genome_End[i] + offset)
    
    if (start_idx > end_idx) {
      print(paste("Skipping ORF", i, "due to invalid index"))
      next
    }
    
    # Extract ORF data
    orf_data <- matrix_data[, start_idx:end_idx]
    
    # Translate to AA
    AA_matrices[[i]] <- t(apply(orf_data, 1, translate))
    
    # Ensure ORF length is a multiple of 3
    orf_length <- end_idx - start_idx + 1
    if (orf_length %% 3 != 0) {
      print(paste("Truncating ORF", i, "to nearest multiple of 3"))
      end_idx <- start_idx + (orf_length %/% 3) * 3 - 1
      orf_data <- matrix_data[, start_idx:end_idx]
    }
    
    # Compute codon matrix size
    num_codons <- (end_idx - start_idx + 1) / 3
    reshaped_orf <- array(orf_data, dim = c(nrow(matrix_data), 3, num_codons))
    codon_matrices[[i]] <- apply(reshaped_orf, c(1, 3), function(x) paste0(x, collapse = ""))
    
    # Assign row names
    rownames(AA_matrices[[i]]) <- names(anc_bayes)
    rownames(codon_matrices[[i]]) <- names(anc_bayes)
  }
  
  return(list(AA_matrices = AA_matrices, codon_matrices = codon_matrices))
}


data_orfs_tb <- read_xlsx("1_Data/1_5_DENV/DenMutationsEpiPosStica2022.xlsx", sheet = "DEN Addresses")
data_orfs_filtered_1 <- dplyr::filter(data_orfs_tb, Serotype == "DENV1")
data_orfs_filtered_2 <- dplyr::filter(data_orfs_tb, Serotype == "DENV2")
data_orfs_filtered_3 <- dplyr::filter(data_orfs_tb, Serotype == "DENV3")
data_orfs_filtered_4 <- dplyr::filter(data_orfs_tb, Serotype == "DENV4")


Codon_AA_output_DENV1 <- process_ORFs_into_AA_and_Codons(
  data_orfs = data_orfs_filtered_1,
  matrix_data = matrix_data_1,
  anc_bayes = anc.bayes_1,
  offset = -94
)


names(Codon_AA_output_DENV1$AA_matrices) <-data_orfs_filtered_1$Protein
names(Codon_AA_output_DENV1$codon_matrices) <-data_orfs_filtered_1$Protein

# Extract the lists from the returned output
list_matrices_ORFs_AA_1 <- Codon_AA_output_DENV1$AA_matrices
list_matrices_ORFs_codon_1 <- Codon_AA_output_DENV1$codon_matrices


## Codon Bounds
# Use codon matrix list (e.g. for DENV1)
codon_list <- list_matrices_ORFs_codon_4

codon_bounds <- lapply(codon_list, function(mat) {
  if (!is.null(mat) && ncol(mat) >= 1) {
    first_codon <- mat[1, 1]
    last_codon <- mat[1, ncol(mat)]
    return(data.frame(First = first_codon, Last = last_codon))
  } else {
    return(data.frame(First = NA, Last = NA))
  }
})

# Combine and format
codon_bounds_df <- do.call(rbind, codon_bounds)
codon_bounds_df$ORF <- names(codon_bounds)
rownames(codon_bounds_df) <- NULL
codon_bounds_df <- codon_bounds_df[, c("ORF", "First", "Last")]


codon_bounds_df


### DENV2-4
Codon_AA_output_DENV2 <- process_ORFs_into_AA_and_Codons(
  data_orfs = data_orfs_filtered_2,
  matrix_data = matrix_data_2,
  anc_bayes = anc.bayes_2,
  offset = -96
)

names(Codon_AA_output_DENV2$AA_matrices) <-data_orfs_filtered_2$Protein
names(Codon_AA_output_DENV2$codon_matrices) <-data_orfs_filtered_2$Protein

# Extract the lists from the returned output
list_matrices_ORFs_AA_2 <- Codon_AA_output_DENV2$AA_matrices
list_matrices_ORFs_codon_2 <- Codon_AA_output_DENV2$codon_matrices




Codon_AA_output_DENV3 <- process_ORFs_into_AA_and_Codons(
  data_orfs = data_orfs_filtered_3,
  matrix_data = matrix_data_3,
  anc_bayes = anc.bayes_3,
  offset = -94
)

names(Codon_AA_output_DENV3$AA_matrices) <-data_orfs_filtered_3$Protein
names(Codon_AA_output_DENV3$codon_matrices) <-data_orfs_filtered_3$Protein

# Extract the lists from the returned output
list_matrices_ORFs_AA_3 <- Codon_AA_output_DENV3$AA_matrices
list_matrices_ORFs_codon_3 <- Codon_AA_output_DENV3$codon_matrices






Codon_AA_output_DENV4 <- process_ORFs_into_AA_and_Codons(
  data_orfs = data_orfs_filtered_4,
  matrix_data = matrix_data_4,
  anc_bayes = anc.bayes_4,
  offset = -179
)

names(Codon_AA_output_DENV4$AA_matrices) <-data_orfs_filtered_4$Protein
names(Codon_AA_output_DENV4$codon_matrices) <-data_orfs_filtered_4$Protein

# Extract the lists from the returned output
list_matrices_ORFs_AA_4 <- Codon_AA_output_DENV4$AA_matrices
list_matrices_ORFs_codon_4 <- Codon_AA_output_DENV4$codon_matrices

View(list_matrices_ORFs_codon_4$whole_genome)



########################################################################################################################################

########################################################################################################################################
## Build dataframes per orf
########################################################################################################################################
codons_possible_with_n = gtools::permutations(n = 5, r = 3, v = c('a', 'c', 't', 'g', 'n'), repeats.allowed = T)
codons_possible_with_n = apply(codons_possible_with_n, MAR = 1, function(x)paste0(x, collapse = ''))
bin_condons <- function(x){
  x=as.numeric(factor(as.character(x),levels=codons_possible_with_n))
  return(x)
}

build_dataframes_node_AA_codons = function(dataset_with_nodes, list_matrices_ORFs_AA, list_matrices_ORFs_codon, orf){
  #' Build Dataframes for Node-Level Amino Acid and Codon Reconstruction
  #'
  #' @description This function processes node-level sequence reconstructions for a given ORF
  #' by extracting amino acid (AA) and codon data from provided ORF matrices. It identifies
  #' single nucleotide polymorphisms (SNPs) in codons and integrates them into the dataset.
  #'
  #' @param dataset_with_nodes Data frame. A dataset containing metadata and node-level sequence reconstructions.
  #' The first four columns should contain metadata, while the remaining data will be matched with ORF matrices.
  #' @param list_matrices_ORFs_AA List. A list of amino acid matrices for different ORFs. Each matrix should have rows
  #' corresponding to sequences and columns representing amino acid positions.
  #' @param list_matrices_ORFs_codon List. A list of codon matrices for different ORFs. Each matrix should have rows
  #' corresponding to sequences and columns representing codon positions.
  #' @param orf Integer or character. The index or name of the ORF to be processed from the input lists.
  #'
  #' @return A list containing two data frames:
  #' \itemize{
  #'   \item{\code{dataset_with_inferred_reconstruction_codon}}{A data frame with metadata and codon SNPs for the selected ORF.}
  #'   \item{\code{dataset_with_inferred_reconstruction_AA}}{A data frame with metadata and amino acid SNPs for the selected ORF.}
  #' }
  #'
  #' @details The function extracts the amino acid and codon sequence data for the specified ORF.
  #' It identifies codon-level SNPs by comparing each sequence to the first sequence in the matrix.
  #' These SNPs are added to the original dataset to create enriched data frames.
  #'
  #' @note This function assumes that sequence names in `dataset_with_nodes` can be extracted from
  #' a path-like string format. If sequence naming conventions differ, modifications may be needed.
  #'
  #' @section SNP Detection:
  #' \subsection{Codon-Level SNPs}{
  #' Codon SNPs are identified by computing the absolute difference between each sequence and
  #' the reference sequence (first row of the codon matrix). Positions with variations are included
  #' in the final dataset.
  #' }
  #'
  #' @examples
  #' # Example Usage:
  #' result <- build_dataframes_node_AA_codons(
  #'   dataset_with_nodes = my_dataset,
  #'   list_matrices_ORFs_AA = list_AA_matrices,
  #'   list_matrices_ORFs_codon = list_codon_matrices,
  #'   orf = 1
  #' )
  #' 
  #' # Extract results
  #' dataset_codon <- result$dataset_with_inferred_reconstruction_codon
  #' dataset_AA <- result$dataset_with_inferred_reconstruction_AA
  #' 
  #' @export
  #' @importFrom stringr str_split
  #'

  # Meta-data with nodes 
  dataset_with_inferred_reconstruction_AA = dataset_with_nodes[,1:4]
  dataset_with_inferred_reconstruction_codon = dataset_with_nodes[,1:4]
  
  # Data AA and codons for that orf
  matrix_ORFs_AA = list_matrices_ORFs_AA[[orf]]
  matrix_ORFs_codon = list_matrices_ORFs_codon[[orf]]
  
  ## Find codon SNPs
  matrix_ORFs_codon_bin = apply(matrix_ORFs_codon, MARGIN = 2, FUN=bin_condons)
  matrix_data_bin_snp = t(apply(matrix_ORFs_codon_bin, MARGIN = 1, FUN=function(x)x-matrix_ORFs_codon_bin[1,]))
  matrix_data_bin_snp = abs(matrix_data_bin_snp)
  overall_snps = colSums(matrix_data_bin_snp)
  a = which(overall_snps > 0)
  
  ## This part is need to align the names of the matrix to the names in dataset_with_inferred_reconstruction_names
  ## YOU MUST CHANGE THIS PART TO ALIGN WITH YOUR NAMING CONVENTION 
  
  # Logical vector: is this entry numeric-looking?
  is_numeric <- grepl("^\\d+$", dataset_with_inferred_reconstruction_codon$name_seq)
  # For numeric entries (assumed to be row indices)
  node_indices <- as.numeric(dataset_with_inferred_reconstruction_codon$name_seq[is_numeric])
  node_names <- rownames(matrix_ORFs_codon)[node_indices]
  # For non-numeric entries (assumed to be tip labels)
  tip_names <- dataset_with_inferred_reconstruction_codon$name_seq[!is_numeric]
  # Merge both into a unified name vector
  dataset_with_inferred_reconstruction_names <- character(nrow(dataset_with_inferred_reconstruction_codon))
  # Currently the outlier is "1867", which doesn’t correspond to any row name in matrix_ORFs_codon.
  # The reconstruction matrix includes 1867 rows.
  # matrix_ORFs_codon only has 1866 rows.
  # "1867" is probably just the final internal node (or an extra dummy tip) not represented in your original codon matrix.
  valid_idx <- idx[!is.na(idx)]
  dataset_with_inferred_reconstruction_names[is_numeric] <- node_names
  dataset_with_inferred_reconstruction_names[!is_numeric] <- tip_names
  
  idx <- match(dataset_with_inferred_reconstruction_names, rownames(matrix_ORFs_codon))
  summary(idx)
  
  
  dataset_with_inferred_reconstruction_codon = cbind(dataset_with_inferred_reconstruction_codon, matrix_ORFs_codon[idx,a])
  dataset_with_inferred_reconstruction_AA = cbind(dataset_with_inferred_reconstruction_AA, matrix_ORFs_AA[idx,a])
  
  colnames(dataset_with_inferred_reconstruction_codon) = c(colnames(dataset_with_inferred_reconstruction_codon)[1:4], a)
  colnames(dataset_with_inferred_reconstruction_AA) = c(colnames(dataset_with_inferred_reconstruction_AA)[1:4], a)
  
  return(list('dataset_with_inferred_reconstruction_codon' = dataset_with_inferred_reconstruction_codon,
              'dataset_with_inferred_reconstruction_AA' = dataset_with_inferred_reconstruction_AA))
  
  
}





## Go through all orfs
for(orf in 1:length(list_matrices_ORFs_AA_1)){
  print(orf)
  reconstruction_sig_tmp = build_dataframes_node_AA_codons(dataset_with_nodes = dataset_with_nodes_1, 
                                                           list_matrices_ORFs_AA = list_matrices_ORFs_AA_1,
                                                           list_matrices_ORFs_codon = list_matrices_ORFs_codon_1,
                                                           orf = orf)
  saveRDS(reconstruction_sig_tmp, file = paste0("4_Output_Data/4_5_DENV/DENV1_Thai_timetree_reconstruction_", data_orfs_filtered_1$Protein[orf], "04_02_25", '.rds'))
}



## Go through all orfs for DENV 2-4
for(orf in 1:length(list_matrices_ORFs_AA_2)){
  print(orf)
  reconstruction_sig_tmp = build_dataframes_node_AA_codons(dataset_with_nodes = dataset_with_nodes_2, 
                                                           list_matrices_ORFs_AA = list_matrices_ORFs_AA_2,
                                                           list_matrices_ORFs_codon = list_matrices_ORFs_codon_2,
                                                           orf = orf)
  saveRDS(reconstruction_sig_tmp, file = paste0("4_Output_Data/4_5_DENV/DENV2_Thai_timetree_reconstruction_", data_orfs_filtered_2$Protein[orf], "04_02_25",  '.rds'))
}

for(orf in 1:length(list_matrices_ORFs_AA_3)){
  print(orf)
  reconstruction_sig_tmp = build_dataframes_node_AA_codons(dataset_with_nodes = dataset_with_nodes_3, 
                                                           list_matrices_ORFs_AA = list_matrices_ORFs_AA_3,
                                                           list_matrices_ORFs_codon = list_matrices_ORFs_codon_3,
                                                           orf = orf)
  saveRDS(reconstruction_sig_tmp, file = paste0("4_Output_Data/4_5_DENV/DENV3_Thai_timetree_reconstruction_", data_orfs_filtered_3$Protein[orf], "04_02_25", '.rds'))
}

for(orf in 1:length(list_matrices_ORFs_AA_4)){
  print(orf)
  reconstruction_sig_tmp = build_dataframes_node_AA_codons(dataset_with_nodes = dataset_with_nodes_4, 
                                                           list_matrices_ORFs_AA = list_matrices_ORFs_AA_4,
                                                           list_matrices_ORFs_codon = list_matrices_ORFs_codon_4,
                                                           orf = orf)
  saveRDS(reconstruction_sig_tmp, file = paste0("4_Output_Data/4_5_DENV/DENV4_Thai_timetree_reconstruction_", data_orfs_filtered_4$Protein[orf], "04_02_25", '.rds'))
}


########################################################################################################################################

# save.image(file = "6_Session_Data/2_5_reconstruct_AA_changes_trees_v2_DAF_04012025img.Rdata")


