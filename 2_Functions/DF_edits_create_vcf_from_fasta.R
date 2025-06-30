
library(vDiveR)
library(tidyr)
library(dplyr)
library(stringr)


# # Example: define gene start-end ranges (in amino acid positions)
# gene_positions <- list(
#   C = c(start = 1, end = 342),
#   prM = c(start = 343, end = 840),
#   E = c(start = 841, end = 2325),
#   NS1 = c(start = 2326, end = 3381),
#   NS2A = c(start = 3382, end = 4038), #Some insertions
#   NS2B = c(start = 4039, end = 4428),
#   NS3 = c(start = 4429, end = 6285),
#   NS4A = c(start = 6286, end = 6666),
#   frag2k = c(start = 6667, end = 6735),
#   NS4B = c(start = 6736, end = 7482),
#   NS5 = c(start = 7483, end = 10188),
#   whole_genome = c(start = 1, end = 10188)
# )
# 
# file_mapping <- list(
#   "1" = list(
#     fasta = here("1_Data", "1_5_DENV", "Seq_vcf_data", "denv_seqs_wgs", "whole_genome", "d1_n1026_underscore.fas"),
#     csv   = here("1_Data", "1_5_DENV", "Seq_vcf_data", "denv_seqs_wgs", "whole_genome", "d1_n1026_geneious_metadata.csv")
#   ),
#   "2" = list(
#     fasta = here("1_Data", "1_5_DENV", "Seq_vcf_data", "denv_seqs_wgs", "whole_genome", "d2_n796_underscore.fas"),
#     csv   = here("1_Data", "1_5_DENV", "Seq_vcf_data", "denv_seqs_wgs", "whole_genome", "d2_n796_geneious_metadata.csv")
#   ),
#   "3" = list(
#     fasta = here("1_Data", "1_5_DENV", "Seq_vcf_data", "denv_seqs_wgs", "whole_genome", "d3_n625_degap_underscore.fas"),
#     csv   = here("1_Data", "1_5_DENV", "Seq_vcf_data", "denv_seqs_wgs", "whole_genome", "d3_n625_degap_geneious_metadata.csv")
#   ),
#   "4" = list(
#     fasta = here("1_Data", "1_5_DENV", "Seq_vcf_data", "denv_seqs_wgs", "whole_genome", "d4_n477_underscore.fas"),
#     csv   = here("1_Data", "1_5_DENV", "Seq_vcf_data", "denv_seqs_wgs", "whole_genome", "d4_n477_geneious_metadata.csv")
#   )
# )
# 
# for (serotype in 1:4) {
#   fasta_path <- file_mapping[[as.character(serotype)]]$fasta
#   out_vcf_dir  <- here("1_Data", "1_5_DENV", "Seq_vcf_data", "denv_seqs_wgs", "vcf_from_fasta", "Test")
#   
#   create_vcf_from_fasta(fasta_path = fasta_path, 
#                         gene_positions = gene_positions, 
#                         output_path = out_vcf_dir, 
#                         serotype = serotype)
#   
# }


create_vcf_from_fasta <- function(fasta_path, gene_positions, output_path, serotype) {
  
  
cat("\n====================\nRunning Serotype", serotype, "\n====================\n")
  



fasta_data = read.fasta(fasta_path)
data_seq = getSequence(fasta_data) #Extracts sequences from the FASTA file.
data_time_seq = str_split_fixed(names(fasta_data), "_", 4)[, 1] #Extract names, which are actually the times
annotations = getAnnot(fasta_data) #Extract the seq ID stored in annotations
data_name_seq <- gsub("^>", "", annotations) # remove the > 


meta_data <- data.frame(
  ID = data_name_seq,
  time = as.numeric(data_time_seq),
  sequence = sapply(data_seq, paste0, collapse = ""),
  Region = str_split_fixed(data_name_seq, "_", 4)[, 2],
  Serotype = str_extract(str_split_fixed(data_name_seq, "_", 4)[, 3], "^\\d+"),
  Clade = str_split_fixed(data_name_seq, "_", 4)[, 3],
  Sample_ID = str_split_fixed(data_name_seq, "_", 4)[, 4]
)


## Transform into matrix of characters
## each row represents a sequence, and each column represents a nucleotide position.
matrix_data = matrix(0, length(data_name_seq), length(data_seq[[1]]))
for (i in 1:length(data_seq)){
  matrix_data[i,] =  data_seq[[i]]
}
rownames(matrix_data) = data_name_seq
remove(fasta_data)

## Put the oldest sequence in the first row: to compute SNP based on it
a = meta_data[which.min(meta_data$time),]
matrix_data = rbind(matrix_data[which(rownames(matrix_data) == a$ID),],
                    matrix_data[-which(rownames(matrix_data) == a$ID),])
rownames(matrix_data) = c(a$ID, rownames(matrix_data)[-1])


# Split matrix_data into a list of matrices per gene
gene_matrices <- lapply(gene_positions, function(range) {
  matrix_data[, range["start"]:range["end"], drop = FALSE]
})



## Transform matrix into binary matrix
bin <- function(x){
  x=as.numeric(factor(as.character(x),levels=c("a","g","c","t", "n", "-", "b", "d", "h", "v", "r", "y", "s", "w", "k", "m")))
  return(x)
}
gene_matrices_bin <- lapply(gene_matrices, function(mat) {
  apply(mat, MARGIN = 2, FUN = bin)
})

for (gene_name in names(gene_matrices_bin)) {
  message("Processing gene: ", gene_name)
  
  matrix_data_bin <- gene_matrices_bin[[gene_name]]
  
  ## Step 1: Compute SNP matrix
  matrix_data_bin_snp <- t(apply(matrix_data_bin, 1, function(x) x - matrix_data_bin[1, ]))
  matrix_data_bin_snp <- apply(matrix_data_bin_snp, 2, function(x) {
    as.numeric(factor(x, levels = c(x[1], -20:(x[1]-1), (x[1]+1):20))) - 1
  })
  
  ## Step 2: Subset to columns with variation
  overall_snps <- colSums(matrix_data_bin_snp)
  snp_cols <- which(overall_snps > 0)
  pos_snp <- snp_cols + gene_positions[[gene_name]]["start"] - 1
  matrix_data_bin_snp_only <- matrix_data_bin_snp[, snp_cols, drop = FALSE]
  
  ## Step 3: Extract single SNPs
  matrix_data_bin_snp_only_singles <- NULL
  pos_snp_singles <- NULL
  pos_snp_singles_names <- NULL
  
  for (i in seq_len(ncol(matrix_data_bin_snp_only))) {
    t <- table(matrix_data_bin_snp_only[, i])
    if (length(t) <= 2) {
      matrix_data_bin_snp_only_singles <- if (is.null(matrix_data_bin_snp_only_singles)) {
        matrix(matrix_data_bin_snp_only[, i], ncol = 1)
      } else {
        cbind(matrix_data_bin_snp_only_singles, matrix_data_bin_snp_only[, i])
      }
      pos_snp_singles <- c(pos_snp_singles, pos_snp[i])
      pos_snp_singles_names <- c(pos_snp_singles_names, pos_snp[i])
    } else {
      for (j in 2:length(t)) {
        tmp <- matrix_data_bin_snp_only[, i]
        tmp[tmp != (j - 1)] <- 0
        tmp[tmp > 0] <- 1
        matrix_data_bin_snp_only_singles <- cbind(matrix_data_bin_snp_only_singles, tmp)
      }
      pos_snp_singles_names <- c(
        pos_snp_singles_names,
        paste0(rep(pos_snp[i], length(t) - 1), "_", LETTERS[1:(length(t) - 1)])
      )
      pos_snp_singles <- c(pos_snp_singles, rep(pos_snp[i], length(t) - 1))
    }
  }
  
  if (is.null(matrix_data_bin_snp_only_singles)) {
    warning("No SNPs found in gene: ", gene_name)
    next
  }
  
  colnames(matrix_data_bin_snp_only_singles) <- pos_snp_singles
  
  ## Step 4: Build VCF
  vcf_header <- c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
  )
  
  data_vcf <- data.frame(
    CHROM = rep(gene_name, ncol(matrix_data_bin_snp_only_singles)),
    POS = pos_snp_singles_names,
    ID = rep(".", ncol(matrix_data_bin_snp_only_singles)),
    REF = rep("NA", ncol(matrix_data_bin_snp_only_singles)),
    ALT = rep("NA", ncol(matrix_data_bin_snp_only_singles)),
    QUAL = rep("5000", ncol(matrix_data_bin_snp_only_singles)),
    FILTER = rep("PASS", ncol(matrix_data_bin_snp_only_singles)),
    INFO = rep(".", ncol(matrix_data_bin_snp_only_singles)),
    FORMAT = rep("GT", ncol(matrix_data_bin_snp_only_singles))
  )
  
  genotype_matrix <- t(matrix_data_bin_snp_only_singles)
  colnames(genotype_matrix) <- data_name_seq
  data_vcf <- cbind(data_vcf, genotype_matrix)
  
  ## Step 5: Write VCF
  gene_output_path <- file.path(output_path, paste0("serotype_", serotype, "_", gene_name, ".vcf"))
  writeLines(vcf_header, con = gene_output_path)
  cat(paste0("#", paste(colnames(data_vcf), collapse = "\t")), file = gene_output_path, append = TRUE, sep = "\n")
  write.table(data_vcf, file = gene_output_path, append = TRUE, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  
  message("VCF written: ", gene_output_path)
}

cat("\n====================\n VCF Output complete for Serotype", serotype, "\n====================\n")
}








