
library(vDiveR)
library(tidyr)
library(dplyr)
library(stringr)


#  Example Usage

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
#   meta_path <- file_mapping[[as.character(serotype)]]$csv
#   out_vcf_path  <- here("1_Data", "1_5_DENV", "Seq_vcf_data", "denv_seqs_wgs", "vcf_from_fasta", paste0("d", serotype, ".vcf"))
#   
#   create_vcf_from_fasta(fasta_path = fasta_path, csv_path = meta_path, output_path = out_vcf_path, serotype = serotype)
# }




create_vcf_from_fasta <- function(fasta_path, csv_path, output_path, serotype) {
  
cat("\n====================\nRunning Serotype", serotype, "\n====================\n")
  
  


meta_data <- read.csv(csv_path)
colnames(meta_data)[colnames(meta_data) == "Name"] <- "time"
colnames(meta_data)[colnames(meta_data) == "Description"] <- "ID"
meta_data_time <-meta_data$time
meta_data <- meta_data %>%
  mutate(
    Region = str_split_fixed(ID, " ", 3)[, 1],
    Clade = str_split_fixed(ID, " ", 3)[, 2],
    Sample_ID = str_split_fixed(ID, " ", 3)[, 3]
  )
meta_data <- meta_data %>% unite(col = "ID", time, ID, sep = " ") %>%
  mutate(
    ID = gsub(" ", "_", ID),
    time = meta_data_time
  )

fasta_data = read.fasta(fasta_path)
data_seq = getSequence(fasta_data) #Extracts sequences from the FASTA file.
data_time_seq = names(fasta_data) #Extract names, which are actually the times
annotations = getAnnot(fasta_data) #Extract the seq ID stored in annotations
data_name_seq <- gsub("^>", "", annotations) # remove the > 


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

## Transform matrix into binary matrix
bin <- function(x){
  x=as.numeric(factor(as.character(x),levels=c("a","g","c","t", "n", "-", "b", "d", "h", "v", "r", "y", "s", "w", "k", "m")))
  return(x)
}
matrix_data_bin = apply(matrix_data, MARGIN = 2, FUN=bin)

## Find SNPs Computes SNPs by comparing each sequence to the reference (first row).
matrix_data_bin_snp = t(apply(matrix_data_bin, MARGIN = 1, FUN=function(x)x-matrix_data_bin[1,]))
matrix_data_bin_snp = apply(matrix_data_bin_snp, MARGIN = 2, FUN=function(x)as.numeric(factor(as.numeric(factor(x, levels = c(x[1], -20:(x[1]-1), (x[1]+1):(20))))))-1)
matrix_data_bin_snp[1:10, 1:10]

## Subset matrix: Retains only positions where at least one SNP is present
overall_snps = colSums(matrix_data_bin_snp)
a = which(overall_snps > 0)
pos_snp = (1:length(data_seq[[1]]))[a]
matrix_data_bin_snp_only = matrix_data_bin_snp[,a]
matrix_data_bin_snp_only[1:10, 1:10]
dim(matrix_data_bin_snp_only)

## Find positions where there are 1+ SNPs
matrix_data_bin_snp_only_singles = NULL
pos_snp_singles = NULL
pos_snp_singles_names = NULL

#Iterates over SNP positions to process single and multiple alleles
for(i in 1:ncol(matrix_data_bin_snp_only)){
  if(i %% 100 == 0){
    print(paste0(i, ' / ', ncol(matrix_data_bin_snp_only)))
  }
  
  t = table(matrix_data_bin_snp_only[,i])
  if(length(t) <= 2){
    if (is.null(matrix_data_bin_snp_only_singles)) {
      matrix_data_bin_snp_only_singles <- matrix(matrix_data_bin_snp_only[, i], ncol = 1)
    } else {
      matrix_data_bin_snp_only_singles <- cbind(matrix_data_bin_snp_only_singles, matrix_data_bin_snp_only[, i])
    }
    pos_snp_singles = c(pos_snp_singles, pos_snp[i])
    pos_snp_singles_names = c(pos_snp_singles_names, pos_snp[i])
    # print(length(pos_snp_singles))
    # print(dim(matrix_data_bin_snp_only_singles))
  }
  if(length(t) > 2){
    # print(t)
    for(j in 2:length(t)){
      tmp = matrix_data_bin_snp_only[,i]
      tmp[which(tmp != (j-1))] = 0
      tmp[which(tmp > 0)] = 1
      matrix_data_bin_snp_only_singles = cbind(matrix_data_bin_snp_only_singles, tmp)
    }
    pos_snp_singles_names = c(pos_snp_singles_names, paste0(rep(pos_snp[i], length(t)-1), '_', LETTERS[1:length(t)-1]))
    pos_snp_singles = c(pos_snp_singles, rep(pos_snp[i], length(t)-1))
  }
}
colnames(matrix_data_bin_snp_only_singles) = pos_snp_singles


## Write genotype vcf
# Create VCF metadata header (manually)
vcf_header <- c(
  "##fileformat=VCFv4.2",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
)

# Prepare VCF table
data_vcf <- data.frame(
  CHROM = rep("CHR1", ncol(matrix_data_bin_snp_only_singles)),
  POS = pos_snp_singles_names,
  ID = rep(".", ncol(matrix_data_bin_snp_only_singles)),
  REF = rep("NA", ncol(matrix_data_bin_snp_only_singles)),
  ALT = rep("NA", ncol(matrix_data_bin_snp_only_singles)),
  QUAL = rep("5000", ncol(matrix_data_bin_snp_only_singles)),
  FILTER = rep("PASS", ncol(matrix_data_bin_snp_only_singles)),
  INFO = rep(".", ncol(matrix_data_bin_snp_only_singles)),
  FORMAT = rep("GT", ncol(matrix_data_bin_snp_only_singles))
)

# Add genotype columns for each sample (rows of matrix)
genotype_matrix <- t(matrix_data_bin_snp_only_singles)
colnames(genotype_matrix) <- data_name_seq  # sample names as column headers
data_vcf <- cbind(data_vcf, genotype_matrix)


# Write header and data to file
writeLines(vcf_header, con = output_path)  # write VCF meta header
write.table(data_vcf, file = output_path, append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
# for i in seq_jc_uni1000_*; do cat line.txt $i > V_$i; done

cat("\n====================\n VCF Output complete for Serotype", serotype, "\n====================\n")
}





