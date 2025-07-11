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

## World only for now

########################################################################################################################################
## Load data
########################################################################################################################################

load("1_Data/1_1_SARS_CoV_2/Initial_index_computation_and_parameters_11102023.Rdata")
load("1_Data/1_1_SARS_CoV_2/Lineages_detected_28062023.Rdata")
## Input fasta data (PHYDAT format)
tree_sars_world_seqs <- read.phyDat("1_Data/1_1_SARS_CoV_2/nextstrain_ncov_gisaid_global_all-time_timetree_dedup_wuhan.fasta",
                                    format = "fasta")
## Input fasta data (SEQINR format)
fasta_data = read.fasta("1_Data/1_1_SARS_CoV_2/nextstrain_ncov_gisaid_global_all-time_timetree_dedup_wuhan.fasta", forceDNAtolower = T, whole.header = T)
data_seq = getSequence(fasta_data)
data_name_seq = names(fasta_data)
########################################################################################################################################

########################################################################################################################################
## Prepare data: make the IDs match
########################################################################################################################################]
correspondence_ID_EPI_world = read.csv("1_Data/1_1_SARS_CoV_2/List_all_isolates_SARSCoV2_global_20230414.csv", 
                                       header = T, col.names = c("N", "Virus name", "Accession ID", "Collection Time", "Country", "Nextstrain_clade", "Source"))
tree = tree_sars_world
names_seqs_world_simple = unlist(lapply(tree$tip.label, function(x){
  tmp = str_split(x, pattern = '\\.')[[1]]
  if(length(tmp) == 1){
    return(tmp)
  }else if(length(tmp) == 2){
    return(tmp[1])
  }else {
    tmp = tmp[-length(tmp)]
    return(paste0(tmp, collapse = '.'))
  }
}))
names(tree_sars_world_seqs) = tree$tip.label[match(correspondence_ID_EPI_world$Virus.name[match(names(tree_sars_world_seqs), correspondence_ID_EPI_world$Accession.ID)], names_seqs_world_simple)]

# tree = tree_sars_world
########################################################################################################################################

########################################################################################################################################
## Reconstruct ancestral sequences within tree
########################################################################################################################################
fit <- pml(tree = tree, data = tree_sars_world_seqs)
fit <- optim.pml(fit, model="GTR", control = pml.control(trace=0))

# makeNodeLabel(tree)

anc.bayes <- ancestral.pml(fit, "bayes")
########################################################################################################################################

########################################################################################################################################
## Write down matrix of sequences (including nodes)
########################################################################################################################################
matrix_data_bin = matrix_data = matrix(0, length(names(anc.bayes)), length(data_seq[[1]]))
for(i in 1:length(names(anc.bayes))){
  matrix_data_bin[i,] =  apply(anc.bayes[[i]], MAR = 1, function(x){
    x = x/sum(x)
    tmp = max(x)
    tmp_which = which.max(x)
    if(tmp > 0.9){return(tmp_which)}else{return(NA)}}
  )[attr(anc.bayes, "index")]
  tmp = factor(matrix_data_bin[i,], levels = 1:4)
  levels(tmp) = attr(anc.bayes, 'levels')
  matrix_data[i,] = as.character(tmp)
  matrix_data[i,which(is.na(matrix_data[i,]))] = 'n'
}
rownames(matrix_data) = names(anc.bayes)

remove(tree_sars_world_seqs) 
remove(fasta_data) 
remove(data_seq)
########################################################################################################################################

########################################################################################################################################
## Check for 1 position that everything is fine
########################################################################################################################################

plot.phylo(tree_sars_world, show.tip.label = F, node.color = matrix_data_bin[,25000])
########################################################################################################################################

########################################################################################################################################
## Split this matrix into ORFs, split by codons and translate to AA
########################################################################################################################################
# data_orfs = read.csv("1_Data/1_1_SARS_CoV_2/Position_ORFs_nextstrain.csv")
# 
# list_matrices_ORFs_AA = list_matrices_ORFs_codon = list()
# for(i in 1:nrow(data_orfs)){
#   if(i==6) {next }
#   i =7 
#   ## AA sequences
#   list_matrices_ORFs_AA[[i]] = t(apply(matrix_data[,data_orfs$start[i]:data_orfs$end[i]], MARGIN = 1, function(x)translate(x)))
#   ## Matrix of codons
#   list_matrices_ORFs_codon[[i]] = matrix(NA, nrow = nrow(matrix_data),
#                                          ncol = ncol(matrix_data[,data_orfs$start[i]:data_orfs$end[i]])/3)
#   for(j in 1:(ncol(list_matrices_ORFs_codon[[i]]))){
#     list_matrices_ORFs_codon[[i]][,j] = apply(matrix_data[,data_orfs$start[i]:data_orfs$end[i]][,(1:3)+(j-1)*(3)], 
#                                               MARGIN = 1, function(x)paste0(x, collapse = ''))
#   }
#   rownames(list_matrices_ORFs_AA[[i]]) = names(anc.bayes)
#   rownames(list_matrices_ORFs_codon[[i]]) = names(anc.bayes)
# }

### Optimized code
data_orfs = read.csv("1_Data/1_1_SARS_CoV_2/Position_ORFs_nextstrain.csv")

# Pre-allocate lists
list_matrices_ORFs_AA = vector("list", length = nrow(data_orfs))
list_matrices_ORFs_codon = vector("list", length = nrow(data_orfs))

# Iterate over ORFs
for (i in seq_len(nrow(data_orfs))) {
  print(paste("Processing ORF", i, "of", nrow(data_orfs)))  # Track progress
  start_idx <- data_orfs$start[i]
  end_idx <- data_orfs$end[i]
  
  # Extract ORF data once to avoid repeated indexing
  orf_data <- matrix_data[, start_idx:end_idx]
  
  # AA sequences (optimized apply usage)
  list_matrices_ORFs_AA[[i]] <- t(apply(orf_data, 1, translate))
  
  # Precompute codon matrix size
  num_codons <- (end_idx - start_idx + 1) / 3
  list_matrices_ORFs_codon[[i]] <- matrix(NA, nrow = nrow(matrix_data), ncol = num_codons)
  
  # Efficient codon extraction using matrix slicing
  for (j in seq_len(num_codons)) {
    col_indices <- (1:3) + (j - 1) * 3
    list_matrices_ORFs_codon[[i]][, j] <- apply(orf_data[, col_indices], 1, paste0, collapse = "")
  }
  
  # Assign row names once (avoid redundant operations)
  rownames(list_matrices_ORFs_AA[[i]]) <- names(anc.bayes)
  rownames(list_matrices_ORFs_codon[[i]]) <- names(anc.bayes)
}


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
  # Meta-data with nodes 
  dataset_with_inferred_reconstruction_AA = dataset_with_nodes[,1:4]
  dataset_with_inferred_reconstruction_codon = dataset_with_nodes[,1:4]
  
  # Data AA and codons fot that orf
  matrix_ORFs_AA = list_matrices_ORFs_AA[[orf]]
  matrix_ORFs_codon = list_matrices_ORFs_codon[[orf]]
  
  ## Find codon SNPs
  matrix_ORFs_codon_bin = apply(matrix_ORFs_codon, MARGIN = 2, FUN=bin_condons)
  matrix_data_bin_snp = t(apply(matrix_ORFs_codon_bin, MARGIN = 1, FUN=function(x)x-matrix_ORFs_codon_bin[1,]))
  matrix_data_bin_snp = abs(matrix_data_bin_snp)
  overall_snps = colSums(matrix_data_bin_snp)
  a = which(overall_snps > 0)
  
  dataset_with_inferred_reconstruction_names = unlist(lapply(dataset_with_inferred_reconstruction_codon$name_seq, function(x)str_split(x, '\\/')[[1]][length(str_split(x, '\\/')[[1]])]))
  idx = match(dataset_with_inferred_reconstruction_names, rownames(matrix_ORFs_codon))
  
  dataset_with_inferred_reconstruction_codon = cbind(dataset_with_inferred_reconstruction_codon, matrix_ORFs_codon[idx,a])
  dataset_with_inferred_reconstruction_AA = cbind(dataset_with_inferred_reconstruction_AA, matrix_ORFs_AA[idx,a])
  
  colnames(dataset_with_inferred_reconstruction_codon) = c(colnames(dataset_with_inferred_reconstruction_codon)[1:4], a)
  colnames(dataset_with_inferred_reconstruction_AA) = c(colnames(dataset_with_inferred_reconstruction_AA)[1:4], a)
  
  return(list('dataset_with_inferred_reconstruction_codon' = dataset_with_inferred_reconstruction_codon,
              'dataset_with_inferred_reconstruction_AA' = dataset_with_inferred_reconstruction_AA))
}

## Go through all orfs
for(orf in 1:length(list_matrices_ORFs_AA)){
  print(orf)
  reconstruction_sig_tmp = build_dataframes_node_AA_codons(dataset_with_nodes = dataset_with_nodes_world, 
                                                            list_matrices_ORFs_AA = list_matrices_ORFs_AA,
                                                            list_matrices_ORFs_codon = list_matrices_ORFs_codon,
                                                            orf = orf)
  saveRDS(reconstruction_sig_tmp, file = paste0("4_Output_Data/4_1_Sars_CoV/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_", data_orfs$gene[orf], '.rds'))
}
########################################################################################################################################


