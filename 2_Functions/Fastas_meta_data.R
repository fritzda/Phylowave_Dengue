library(seqinr)
library(here)

library(vDiveR)
library(tidyr)
library(dplyr)
library(stringr)



d1_fasta_file= "/Users/fritzda/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/Collaborators_Data/Dengue_Sequences_to_Date_Hat/20250206_aligned/d4.fasta"
d1_fasta_data = read.fasta(d1_fasta_file)
d1_data_seq = getSequence(d1_fasta_data) #Extracts sequences from the FASTA file.
d1_data_time_seq = names(d1_fasta_data) #Extract names, which are actually the times
d1_data_name_seq = getAnnot(d1_fasta_data) #Extract the seq ID stored in annotations
d1_data_name_seq = clean_annotations <- gsub("^>", "", d1_data_name_seq) #Remove the >

## Transform into matrix of characters
## each row represents a sequence, and each column represents a nucleotide position.
d1_matrix_data = matrix(0, length(d1_data_name_seq), length(d1_data_seq[[1]]))
for (i in 1:length(d1_data_seq)){
  d1_matrix_data[i,] =  d1_data_seq[[i]]
}
rownames(d1_matrix_data) = d1_data_name_seq
remove(d1_fasta_data)



meta_data <- as.data.frame(d1_data_name_seq)
meta_data$Long_ID <- meta_data$d1_data_name_seq
meta_data <- meta_data %>%
  mutate(
    Time = str_split_fixed(Long_ID, " ", 4)[,1],
    Region = str_split_fixed(Long_ID, " ", 4)[, 2],
    Clade = str_split_fixed(Long_ID, " ", 4)[, 3],
    Sample_ID = str_split_fixed(Long_ID, " ", 4)[, 4]
  )
meta_data <- meta_data[,2:6]
# meta_data <- meta_data %>% unite(col = "Long", Time, ID, sep = " ")
# meta_data <- meta_data %>%
#   mutate(
#     time = meta_data_time
#   )

outpath <- "/Users/fritzda/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/Collaborators_Data/Dengue_Sequences_to_Date_Hat/20250206_aligned/d4_fasta_DAF_Rmade_metadata.csv"
write.csv(meta_data, outpath)


## Put the oldest sequence in the first row: to compute SNP based on it
a1 = meta_data[which.min(meta_data$Time),]
d1_matrix_data = rbind(d1_matrix_data[which(rownames(d1_matrix_data) == a1$Long_ID),],
                    d1_matrix_data[-which(rownames(d1_matrix_data) == a1$Long_ID),])
rownames(d1_matrix_data) = c(a1$Long_ID, rownames(d1_matrix_data)[-1])

outpath <- "/Users/fritzda/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/Collaborators_Data/Dengue_Sequences_to_Date_Hat/20250206_aligned/d4_fasta_DAF_Rmade_seq_matrix.csv"
write.csv(d1_matrix_data, outpath)
