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



########################################################################################################################################
## Load data
########################################################################################################################################
here::i_am("DENV_3_fitness.Rmd")

# load('2_analysis_index/1_index_computations/Initial_index_computation_and_parameters_11102023.Rdata')
# load('2_analysis_index/2_find_index_groups/Lineages_detected_11102023.Rdata')

## Input fasta data (PHYDAT format)
tree_DENV1 <- read.phyDat("1_Data/1_5_DENV/Seq_vcf_data/denv_seqs_wgs/whole_genome/d1_n1026.fas",
                              format = "fasta")
names(tree_DENV1) <- gsub(" ", "_", names(tree_DENV1))

## Input fasta data (SEQINR format)
fasta_data_DENV1 = read.fasta('1_Data/1_5_DENV/Seq_vcf_data/denv_seqs_wgs/whole_genome/d1_n1026.fas', forceDNAtolower = T, whole.header = T)
data_seq_DENV1 = getSequence(fasta_data_DENV1)
data_name_seq_DENV1 = names(fasta_data_DENV1)
data_name_seq_DENV1 = gsub(" ", "_", data_name_seq_DENV1)



#Read in tree as phylo 
tree_phylo_DENV1 = read.nexus("1_Data/1_5_DENV/best_tree_d1_Lin.nexus")
tree_phylo_DENV1 = collapse.singles(ladderize(multi2di(tree_phylo_DENV1, random = F), right = F))

setdiff( names(tree_DENV1), tree_phylo_DENV1$tip.label)
# 934 tip labels vs 1026 in the treedata


########################################################################################################################################

########################################################################################################################################
## Prepare data: make the IDs match
########################################################################################################################################]
# correspondence_ID_EPI_world = read.csv('1_Data_Nextstrain_20230414/World/nextstrain_ncov_gisaid_global_all-time_timetree_virus_ID.csv', 
#                                        header = F, col.names = c("N", "Virus name", "Accession ID"))
correspondence_ID_EPI_genious_DENV1 = read.csv("1_Data/1_5_DENV/Seq_vcf_data/denv_seqs_wgs/whole_genome/d1_n1026_geneious_metadata.csv", 
                                               col.names = c("Name",	"Created Date",	"Description",	"Imported From: Filename",	"Imported From: Path",	"Sequence", "A", "B"))


# Remove the dates with multiples designation in metadata names  ("2012.361 2" → "2012.361") )
correspondence_ID_EPI_genious_DENV1$Name <- gsub("^([0-9]+\\.[0-9]+).*", "\\1", correspondence_ID_EPI_genious_DENV1$Name) 
  # ^([0-9]+\\.[0-9]+).*
  #  ^ → Start of string.
  #  ([0-9]+\\.[0-9]+) → Captures a pattern of a year followed by a decimal (e.g., "2011.501").
  # .* → Matches anything that follows.
  #  "\\1" → Keeps only the first matched group (the year and decimal part), removing everything else.


correspondence_ID_EPI_all <- list(correspondence_ID_EPI_genious_DENV1)

names(correspondence_ID_EPI_all) <- c("DENV1")

# Apply transformations to all data frames in the list, create longname and make it look liek the vcf file names
correspondence_ID_EPI_all <- lapply(correspondence_ID_EPI_all, function(df) {
  # Create the 'longname' column dynamically
  df$longname <- gsub(" ", "_", paste(df$Name, df$Description, sep = "_"))
  
  # Reorder columns to put 'longname' first
  df <- df[, c("longname", setdiff(names(df), "longname"))]
  
  return(df)
})


# match(correspondence_ID_EPI_all$DENV1$longname,data_name_seq_DENV1)


# tree = tree_sars_world
# names_seqs_world_simple = unlist(lapply(tree$tip.label, function(x){
#   tmp = str_split(x, pattern = '\\.')[[1]] #Split at the periods
#  # if the result has 1 part return as is
#    if(length(tmp) == 1){ 
#     return(tmp)
#  # if it has 2 parts, return only the first
#   }else if(length(tmp) == 2){
#     return(tmp[1])
#  # if it has any other number of parts: remove the last parts and return the rest as a joined string with periods
#   }else {
#     tmp = tmp[-length(tmp)]
#     return(paste0(tmp, collapse = '.'))
#   }
# }))
# #match on virus name and accesion id 
# names(tree_sars_world_seqs) = tree$tip.label[match(correspondence_ID_EPI_world$Virus.name[match(names(tree_sars_world_seqs), correspondence_ID_EPI_world$Accession.ID)], names_seqs_world_simple)]


########################################################################################################################################

########################################################################################################################################
## Reconstruct ancestral sequences within tree
########################################################################################################################################
tree = tree_phylo_DENV1

# pml() A phangorn package that create a likelyhood of a tree, requires a phylo class tree and then a phydat data object
# optim.pml allows using a GTR model for evolution
fit <- pml(tree = tree, data = tree_DENV1)
fit <- optim.pml(fit, model="GTR", control = pml.control(trace=0))

setdiff( names(tree_DENV1), names(fit$data))



# Marginal reconstruction of the ancestral character states.
anc.bayes <- ancestral.pml(fit, "bayes")


########################################################################################################################################

########################################################################################################################################
## Write down matrix of sequences (including nodes)
########################################################################################################################################

# Creates two matrices, matrix_data_bin and matrix_data, both initialized with zeros.
matrix_data_bin = matrix_data = matrix(0, length(names(anc.bayes)), length(data_seq_DENV1[[1]]))

# #Loop Over Ancestral Sequences
# for(i in 1:length(names(anc.bayes))){ 
#   matrix_data_bin[i,] =  apply(anc.bayes[[i]], MAR = 1, function(x){
#     x = x/sum(x)  # Normalize the probabilities
#     tmp = max(x)  # Find the highest probability
#     tmp_which = which.max(x)  # Find the index of the highest probability
#     if(tmp > 0.9){
#       return(tmp_which)  # If the max probability is greater than 0.9, assign that state
#     }
#     else{
#       return(NA)  # Otherwise, assign NA
#     }
#   }
#   )[attr(anc.bayes, "index")]
#   tmp = factor(matrix_data_bin[i,], levels = 1:4)
#   levels(tmp) = attr(anc.bayes, 'levels')
#   matrix_data[i,] = as.character(tmp)
#   matrix_data[i,which(is.na(matrix_data[i,]))] = 'n'
# }


for(i in 1:length(names(anc.bayes))) {
  matrix_data_bin[i,] = apply(anc.bayes[[i]], MAR = 1, function(x) {
    if (all(is.na(x)) || sum(x, na.rm = TRUE) == 0) {
      return(NA)  # If all values are NA or sum is zero, return NA to avoid errors
    }
    x = x / sum(x, na.rm = TRUE)  # Normalize probabilities safely
    tmp = max(x, na.rm = TRUE)  # Find the max probability safely
    tmp_which = which.max(x)  # Get the index of the most probable nucleotide
    
    if (!is.na(tmp) && tmp > 0.9) {
      return(tmp_which)
    } else {
      return(NA)
    }
  })[attr(anc.bayes, "index")]
  
  tmp = factor(matrix_data_bin[i,], levels = 1:4)
  levels(tmp) = attr(anc.bayes, 'levels')
  matrix_data[i,] = as.character(tmp)
  matrix_data[i, which(is.na(matrix_data[i,]))] = 'n' #Assings n if the value is NA
}

rownames(matrix_data) = names(anc.bayes)

remove(tree_DENV1) 
remove(fasta_data_DENV1) 
remove(data_seq_DENV1)
########################################################################################################################################

########################################################################################################################################
## Check for 1 position that everything is fine, this can take a minute to render
########################################################################################################################################
plot.phylo(tree_phylo_DENV1, show.tip.label = F, node.color = matrix_data_bin[,500])
########################################################################################################################################

########################################################################################################################################
## Split this matrix into ORFs, split by codons and translate to AA
########################################################################################################################################
data_orfs = read.csv('1_refseq/Position_ORFs_nextstrain.csv')

list_matrices_ORFs_AA = list_matrices_ORFs_codon = list()
for(i in 1:nrow(data_orfs)){
  print(i)
  ## AA sequences
  list_matrices_ORFs_AA[[i]] = t(apply(matrix_data[,data_orfs$start[i]:data_orfs$end[i]], MARGIN = 1, function(x)translate(x)))
  ## Matrix of codons
  list_matrices_ORFs_codon[[i]] = matrix(NA, nrow = nrow(matrix_data),
                                         ncol = ncol(matrix_data[,data_orfs$start[i]:data_orfs$end[i]])/3)
  for(j in 1:(ncol(list_matrices_ORFs_codon[[i]]))){
    list_matrices_ORFs_codon[[i]][,j] = apply(matrix_data[,data_orfs$start[i]:data_orfs$end[i]][,(1:3)+(j-1)*(3)], 
                                              MARGIN = 1, function(x)paste0(x, collapse = ''))
  }
  rownames(list_matrices_ORFs_AA[[i]]) = names(anc.bayes)
  rownames(list_matrices_ORFs_codon[[i]]) = names(anc.bayes)
}
########################################################################################################################################

########################################################################################################################################
## Build dataframes per orf
########################################################################################################################################
codons_possibale_with_n = gtools::permutations(n = 5, r = 3, v = c('a', 'c', 't', 'g', 'n'), repeats.allowed = T)
codons_possibale_with_n = apply(codons_possibale_with_n, MAR = 1, function(x)paste0(x, collapse = ''))
bin_condons <- function(x){
  x=as.numeric(factor(as.character(x),levels=codons_possibale_with_n))
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
  saveRDS(reconstruction_sig_tmp, file = paste0('2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_', data_orfs$gene[orf], '.rds'))
}
########################################################################################################################################


