## Packages used
library(stringr)
library(ape)
library(phangorn)
library(BactDating)
library(phytools)
library(coda)
library(thd)
library(vcfR)
library(lubridate)
library(ggplot2)
library(ggtree)

########################################################################################################################################
## Useful functions
########################################################################################################################################
## Plotting functions
mean.and.ci <-function(v){ return( c(mean(v), as.numeric(quantile(v,probs = 0.025, na.rm = T)), as.numeric(quantile(v,probs = 0.975, na.rm = T))))}
axisPhylo_NL = function (side = 1, root.time = NULL, backward = TRUE, at_axis = NULL, lab_axis = NULL, ...){
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  type <- lastPP$type
  if (type == "unrooted")
    stop("axisPhylo() not available for unrooted plots; try add.scale.bar()")
  if (type == "radial")
    stop("axisPhylo() not meaningful for this type of plot")
  if (is.null(root.time))
    root.time <- lastPP$root.time
  if (type %in% c("phylogram", "cladogram")) {
    xscale <- if (lastPP$direction %in% c("rightwards", "leftwards"))
      range(lastPP$xx)
    else range(lastPP$yy)
    tmp <- lastPP$direction %in% c("leftwards", "downwards")
    tscale <- c(0, xscale[2] - xscale[1])
    if (xor(backward, tmp))
      tscale <- tscale[2:1]
    if (!is.null(root.time)) {
      tscale <- tscale + root.time
      if (backward)
        tscale <- tscale - xscale[2]
    }
    beta <- diff(xscale)/diff(tscale)
    alpha <- xscale[1] - beta * tscale[1]
    if(is.null(at_axis) == T){
      x <- beta * lab + alpha
      lab <- pretty(tscale)
    }
    # if(is.null(at_axis) != F){
    x <- at_axis
    lab <- lab_axis
    # }
    axis(side = side, at = x, labels = lab, ...)
  }
  else {
    n <- lastPP$Ntip
    xx <- lastPP$xx[1:n]
    yy <- lastPP$yy[1:n]
    r0 <- max(sqrt(xx^2 + yy^2))
    alpha <- sort(setNames(rect2polar(xx, yy)$angle, 1:n))
    angles <- c(diff(alpha), 2 * pi - alpha[n] + alpha[1L])
    j <- which.max(angles)
    i <- if (j == 1L)
      n
    else j - 1L
    firstandlast <- as.integer(names(angles[c(i, j)]))
    theta0 <- mean(atan2(yy[firstandlast], xx[firstandlast]))
    x0 <- r0 * cos(theta0)
    y0 <- r0 * sin(theta0)
    inc <- diff(pretty(c(0, r0))[1:2])
    srt <- 360 * theta0/(2 * pi)
    coef <- -1
    if (abs(srt) > 90) {
      srt <- srt + 180
      coef <- 1
    }
    len <- 0.025 * r0
    r <- r0
    while (r > 1e-08) {
      x <- r * cos(theta0)
      y <- r * sin(theta0)
      if (len/r < 1) {
        ra <- sqrt(len^2 + r^2)
        thetaa <- theta0 + coef * asin(len/r)
        xa <- ra * cos(thetaa)
        ya <- ra * sin(thetaa)
        segments(xa, ya, x, y)
        text(xa, ya, r0 - r, srt = srt, adj = c(0.5,
                                                1.1), ...)
      }
      r <- r - inc
    }
    segments(x, y, x0, y0)
  }
}
## THD function
thd_home_weights = function(H, t, m, mu, scale = "time", q = 0.5, weights) {
  # Compute THD bandwidth
  b <- thd.bandwidth(t, m, mu, q)
  # Check symmetry
  # stopifnot(all(H == t(H)))
  # Check zero diagonal
  stopifnot(all(diag(H) == 0))
  # Check correction matrix has the same dimensions as H
  stopifnot(all(dim(H) == dim(weights)))
  # Let B <- b^h
  B <- b^H
  # Take columns sums, subtract 1 (= b^0) to ignore diagonal
  # Bsums <- colSums(B, na.rm = T) - 1
  # Take columns sums, of the corrected matrix
  # Bsums_corrected <- colSums(B*weights) ## correction is a matrix of 0s & 1s, same dim as H, 1=sample taken into account, 0=not
  # Take columns sums, of the corrected matrix
  Bsums_corrected2 <- colSums(B, na.rm = T)-1  ## correction is a matrix of 0s & 1s, same dim as H, 1=sample taken into account, 0=not
  # Normalization constant
  # nc <- (1 - b) / ( (1 - b^(m+1)) * (nrow(H) - 1) )
  # Rescale if required
  return(switch(scale,
                density = {
                  # Normalized density
                  Bsums * nc
                },
                size = {
                  # Rescale by no. of markers
                  Bsums * nc * m
                },
                relative = {
                  # Relative to completely homogeneous population
                  Bsums / (nrow(H) - 1)
                },
                corrected = {
                  # Minimizing sampling bias
                  Bsums_corrected/(colSums(weights))
                }, 
                corrected2 = {
                  # Minimizing sampling bias
                  Bsums_corrected2/(colSums(weights)-1)
                }, 
                time = {
                  # Relative to expected distance since timescale
                  Bsums / ( (nrow(H) - 1) * b^iam.distance(t, m, mu))
                }
  ))
}
THD_from_tips_and_nodes_time_windows = function(dist_mat, t, wind, data, mu){
  names_seqs = t$tip.label
  
  ## Construct a new distance matrix, that includes branches (unsampled individuals in the population)
  H_with_sampled = matrix(NA, nrow = nrow(dist_mat), ncol = ncol(dist_mat))
  colnames(H_with_sampled) = colnames(dist_mat)
  rownames(H_with_sampled) = colnames(dist_mat)
  
  ## Compute branches birth and death
  branches_times = t$edge
  branches_times[,1] = branches_times[,2] = NA
  times_tips_nodes = data$time
  a = match(t$edge[,1], data$ID)
  branches_times[,1] = times_tips_nodes[a]
  a = match(t$edge[,2], data$ID)
  branches_times[,2] = times_tips_nodes[a]
  
  ## Loop to correct the matrix for each individual
  for(i in 1:ncol(H_with_sampled)){
    # print(i)
    sample = colnames(H_with_sampled)[i]
    data_tmp = data[which(data$name_seq == sample),]
    
    ## Filter branches alive at sampling time
    t_min = data_tmp$time - wind
    t_max = data_tmp$time + wind
    a = which(branches_times[,1] <= t_min & branches_times[,2] >= t_max) ## branches alive, with no sampled individuals 
    b = which(branches_times[,1] <= t_min & (branches_times[,2] >= t_min & branches_times[,2] <= t_max)) ## branches born before interval and died within interval
    c = which((branches_times[,1] >= t_min & branches_times[,1] <= t_max) & branches_times[,2] >= t_max) ## branches born in interval and died after interval
    d = which((branches_times[,1] >= t_min & branches_times[,1] <= t_max) & (branches_times[,2] >= t_min & branches_times[,2] <= t_max) ) ## branches born and died within interval
    
    ## For each sampling above, get distances
    edge_lengths = branches_times[c(a,b,c,d),]
    edges = t$edge[c(a,b,c,d),]
    
    h_tmp = dist_mat[,i] ## Distances between these chosen tips and nodes
    
    if(is.null(dim(edges))){
      ## List all nodes in the tables
      # all_nodes_and_tips_times = sort(c(edge_lengths[,2])) ## Choice: consider one individual per branch, and the lastest sample (offspring)
      all_nodes_and_tips = sort(c(edges[2])) ## Choice: consider one individual per branch, and the lastest sample (offspring)
      h_tmp[is.na(match(data$ID, unique(all_nodes_and_tips)))] = NA ## Only keep the relevant nodes
      
      ## Self distance should be 0
      h_tmp[i] = 0
    }
    if(!is.null(dim(edges))){
      ## List all nodes in the tables
      # all_nodes_and_tips_times = sort(c(edge_lengths[,2])) ## Choice: consider one individual per branch, and the lastest sample (offspring)
      all_nodes_and_tips = sort(c(edges[,2])) ## Choice: consider one individual per branch, and the lastest sample (offspring)
      h_tmp[is.na(match(data$ID, unique(all_nodes_and_tips)))] = NA ## Only keep the relevant nodes
      
      ## Remove connecting nodes: those are messing up the distance matrix
      tbl = table(c(edges[,1], edges[,2])) 
      ## Any node that is present 3 times should be removed: it's a connecting node - which is screwing up the distance matrix later
      ## Any node that is present 2 times is a mrca, but is already not taken into account is the computation, so it's fine
      connecting_nodes = as.numeric(names(tbl)[which(tbl == 3)])
      m = !is.na(match(data$ID, connecting_nodes))
      h_tmp[m] = NA ## Remove connecting nodes
      
      if(length(which(connecting_nodes == i)) > 0){ ## The node i is a connecting node
        m = match(data$ID, edges[which(edges[,1] == i),2])
        m = which(!is.na(m))
        h_tmp[m] = NA ## Remove tips descending from connecting nodes
      }
      
      ## Correct this vector for the nodes that have a sampling time >data_tmp$time
      too_late = which(edge_lengths[,2] > data_tmp$time)
      h_tmp[match(edges[too_late,2], data$ID)] = h_tmp[match(edges[too_late,2], data$ID)]  - abs(edge_lengths[too_late,2] - data_tmp$time)
      
      ## Correct this vector for the nodes that have a sampling time <data_tmp$time
      too_early = which(edge_lengths[,2] < data_tmp$time)
      # h_tmp[edges[too_early,2]] = h_tmp[edges[too_early,2]]  + abs(edge_lengths[too_early,2] - data_tmp$time)
      h_tmp[match(edges[too_early,2], data$ID)] = h_tmp[match(edges[too_early,2], data$ID)]  + abs(edge_lengths[too_early,2] - data_tmp$time)
      
      ## Self distance should be 0
      h_tmp[i] = 0
    }
    ## Put this vector in the big matrix 
    H_with_sampled[,i] = h_tmp
  }
  
  ## Isolates within the same range
  same_year_mat = apply(H_with_sampled, MARGIN = 2, FUN = function(x)!is.na(x))
  same_year_mat = !is.na(H_with_sampled)
  
  ## Compute THD, by year
  THD = thd_home_weights(H_with_sampled*mu*genome_length/2, timescale, genome_length, mu, scale = 'corrected2', weights = same_year_mat)
  
  return(THD)
}
## Chain extraction
read.chains.from.table = function(table){
  if(typeof(table) != 'double') {
    print('Changing type of table to matrix')
    table = as.matrix(table)
  }
  Chains = list()
  col_variables = sapply(colnames(table), function(x)str_split(x, pattern = "[.]")[[1]][1])
  variable_names = unique(col_variables)
  nchains = nrow(table)
  for(i in 1:length(variable_names)){
    a = match(col_variables, variable_names[i])
    if(length(which(is.na(a) == F)) == 1) {
      Chains[[i]] = table[,match(variable_names[i], col_variables)]
    }
    else{
      a = match(col_variables, variable_names[i])
      tmp = colnames(table)[which(is.na(a) == F)]
      ndims = 0
      dims = NULL
      empty = F
      while(empty == F & ndims < 10){
        l = sapply(tmp, function(x)str_split(x, pattern = "[.]")[[1]][1+ndims+1])
        if(all(is.na(l))) empty = T
        if(all(is.na(l)) == F) {
          ndims = ndims + 1
          dims = c(dims, max(as.numeric(l)))
        }
      }
      if(ndims >5) print('Error: this function only supports arrays of <=5 dimensions')
      Chains[[i]] = array(NA, dim = c(nchains, dims))
      a = which(is.na(match(col_variables, variable_names[i]))==F)
      
      if(ndims == 1){
        Chains[[i]] = table[,a]
        colnames(Chains[[i]]) = NULL
      }
      
      if(ndims == 2){
        j = 1
        for(d2 in 1:dims[2]){
          for(d1 in 1:dims[1]){
            Chains[[i]][,d1,d2] = table[,a[j]]
            j = j+1
          }
        }
      }
      
      if(ndims == 3){
        j = 1
        for(d3 in 1:dims[3]){
          for(d2 in 1:dims[2]){
            for(d1 in 1:dims[1]){
              Chains[[i]][,d1,d2,d3] = table[,a[j]]
              j = j+1
            }
          }
        }
      }
      
      if(ndims == 4){
        j = 1
        for(d4 in 1:dims[4]){
          for(d3 in 1:dims[3]){
            for(d2 in 1:dims[2]){
              for(d1 in 1:dims[1]){
                Chains[[i]][,d1,d2,d3,d4] = table[,a[j]]
                j = j+1
              }
            }
          }
        }
      }
      
      if(ndims == 5){
        j = 1
        for(d5 in 1:dims[5]){
          for(d4 in 1:dims[4]){
            for(d3 in 1:dims[3]){
              for(d2 in 1:dims[2]){
                for(d1 in 1:dims[1]){
                  Chains[[i]][,d1,d2,d3,d4] = table[,a[j]]
                  j = j+1
                }
              }
            }
          }
        }
      }
    }
    names(Chains)[i] = variable_names[i]
  }
  return(Chains)
}
########################################################################################################################################

########################################################################################################################################
## Load data & set parameters
########################################################################################################################################
## Set working directory
setwd('/Users/noemielefrancq/Documents/THD/THD_DENV/')

## Timed-trees
tree_d1 = read.nexus(paste0('Data/d', 1, '.fas.timetree.nexus'))
tree_d1 = multi2di(tree_d1)
tree_d2 = read.nexus(paste0('Data/d', 2, '.fas.timetree.nexus'))
tree_d2 = multi2di(tree_d2)
tree_d3 = read.nexus(paste0('Data/d', 3, '.fas.timetree.nexus'))
tree_d3 = multi2di(tree_d3)
tree_d4 = read.nexus(paste0('Data/d', 4, '.fas.timetree.nexus'))
tree_d4 = multi2di(tree_d4)

## Metadata
meta_data_d1 = read.csv(file = paste0("Data/d", 1, ".csv"), sep = ',', header = T)
meta_data_d2 = read.csv(file = paste0("Data/d", 2, ".csv"), sep = ',', header = T)
meta_data_d3 = read.csv(file = paste0("Data/d", 3, "_degap.csv"), sep = ',', header = T)
meta_data_d4 = read.csv(file = paste0("Data/d", 4, ".csv"), sep = ',', header = T)

## Vcfs nucleotides vcfs
data_vcf_d1 = read.csv(file = paste0("Data/d", 1, '.vcf'), sep = '\t')
data_vcf_d2 = read.csv(file = paste0("Data/d", 2, '.vcf'), sep = '\t')
data_vcf_d3 = read.csv(file = paste0("Data/d", 3, '.vcf'), sep = '\t')
data_vcf_d4 = read.csv(file = paste0("Data/d", 4, '.vcf'), sep = '\t')

## Vcfs prots
data_vcf_anc_cap_C = read.csv(file = 'Data/AA/anchored capsid protein C.vcf', sep = '\t')
data_vcf_cap_C = read.csv(file = 'Data/AA/capsid protein C.vcf', sep = '\t')
data_vcf_E = read.csv(file = 'Data/AA/envelope protein E.vcf', sep = '\t')
data_vcf_gly_M = read.csv(file = 'Data/AA/membrane glycoprotein M.vcf', sep = '\t')
data_vcf_gly_pre_M = read.csv(file = 'Data/AA/membrane glycoprotein precursor M.vcf', sep = '\t')
data_vcf_NS1 = read.csv(file = 'Data/AA/nonstructural protein NS1.vcf', sep = '\t')
data_vcf_NS2A = read.csv(file = 'Data/AA/nonstructural protein NS2A.vcf', sep = '\t')
data_vcf_NS2B = read.csv(file = 'Data/AA/nonstructural protein NS2B.vcf', sep = '\t')
data_vcf_NS3 = read.csv(file = 'Data/AA/nonstructural protein NS3.vcf', sep = '\t')
data_vcf_NS4A = read.csv(file = 'Data/AA/nonstructural protein NS4A.vcf', sep = '\t')
data_vcf_NS4B = read.csv(file = 'Data/AA/nonstructural protein NS4B.vcf', sep = '\t')
data_vcf_prot_2K = read.csv(file = 'Data/AA/protein 2K.vcf', sep = '\t')
data_vcf_prot_pr = read.csv(file = 'Data/AA/protein pr.vcf', sep = '\t')
data_vcf_RNA_pol = read.csv(file = 'Data/AA/RNA-dependent RNA polymerase NS5.vcf', sep = '\t')
## Put all vcfs into list
data_vcf_list = list(data_vcf_anc_cap_C, data_vcf_cap_C, data_vcf_E, data_vcf_gly_M, data_vcf_gly_pre_M,
                     data_vcf_NS1, data_vcf_NS2A, data_vcf_NS2B, data_vcf_NS3, data_vcf_NS4A, data_vcf_NS4B,
                     data_vcf_prot_2K, data_vcf_prot_pr, data_vcf_RNA_pol)
names(data_vcf_list) = c('anc_cap_C', 'cap_C', 'E', 'gly_M', 'gly_pre_M',
                            'NS1', 'NS2A', 'NS2B', 'NS3', 'NS4A', 'NS4B',
                            'prot_2K', 'prot_pr', 'RNA_pol')

## Names all sequences
names_seqs_d1 = tree_d1$tip.label
names_seqs_d2 = tree_d2$tip.label
names_seqs_d3 = tree_d3$tip.label
names_seqs_d4 = tree_d4$tip.label

## Mutation rate from Hat's log files (LSD dating)
mu_d1 = 0.000716792 ## in mutations/site/year
mu_d2 = 0.000428446 ## in mutations/site/year
mu_d3 = 0.000110237 ## in mutations/site/year
mu_d4 = 0.000365477 ## in mutations/site/year

## Length genome 
genome_length = 10188 ## Length genome

## Parameters for THD
timescale = 2 ## Timescale
########################################################################################################################################

########################################################################################################################################
## Preparation data tips
########################################################################################################################################
## Create dataset with names of each sequence, time sampling, prn type etc
add_snp_data = function(dataset_tips, names_seqs, data_vcf_list){
  for(i in 1:length(data_vcf_list)){
    tmp = data_vcf_list[[i]]
    tmp = cbind(tmp, rep(NA, nrow(tmp)))
    a = match(names_seqs, unlist(lapply(colnames(data_vcf_list[[i]]), function(x)strsplit(x, split = 'X')[[1]][2])))
    a[which(is.na(a))] = ncol(tmp)
    df_tmp = t(tmp[,a])
    colnames(df_tmp) = paste0(data_vcf_list[[i]]$REF, data_vcf_list[[i]]$POS, data_vcf_list[[i]]$ALT, '_',names(data_vcf_list)[i])
    dataset_tips = cbind(dataset_tips, df_tmp)
  }
  return(dataset_tips)
}
  
## DENV1
dataset_tips_d1 = data.frame('ID' = 1:length(names_seqs_d1),
                          'name_seq' = names_seqs_d1,
                          'time' = meta_data_d1$time[match(names_seqs_d1, meta_data_d1$shortname)])
dataset_tips_d1 = add_snp_data(dataset_tips = dataset_tips_d1, names_seqs = names_seqs_d1, data_vcf_list)

## DENV2
dataset_tips_d2 = data.frame('ID' = 1:length(names_seqs_d2),
                          'name_seq' = names_seqs_d2,
                          'time' = meta_data_d2$time[match(names_seqs_d2, meta_data_d2$shortname)])
dataset_tips_d2 = add_snp_data(dataset_tips = dataset_tips_d2, names_seqs = names_seqs_d2, data_vcf_list)

## DENV3
dataset_tips_d3 = data.frame('ID' = 1:length(names_seqs_d3),
                          'name_seq' = names_seqs_d3,
                          'time' = meta_data_d3$time[match(names_seqs_d3, meta_data_d3$shortname)])
dataset_tips_d3 = add_snp_data(dataset_tips = dataset_tips_d3, names_seqs = names_seqs_d3, data_vcf_list)

## DENV4
dataset_tips_d4 = data.frame('ID' = 1:length(names_seqs_d4),
                          'name_seq' = names_seqs_d4,
                          'time' = meta_data_d4$time[match(names_seqs_d4, meta_data_d4$shortname)])
dataset_tips_d4 = add_snp_data(dataset_tips = dataset_tips_d4, names_seqs = names_seqs_d4, data_vcf_list)
########################################################################################################################################

########################################################################################################################################
## Compute distance between each pair of sequences and nodes in the tree
########################################################################################################################################
# D1
genetic_distance_mat_d1 = dist.nodes(tree_d1) ## Pairwise time of differences inferred with BEAST (ie taking into account substitution model)
colnames(genetic_distance_mat_d1) = c(names_seqs_d1, length(names_seqs_d1)+(1:(length(names_seqs_d1)-1)))
rownames(genetic_distance_mat_d1) = c(names_seqs_d1, length(names_seqs_d1)+(1:(length(names_seqs_d1)-1)))
hamming_dist_raw_d1 = genetic_distance_mat_d1

# D2
genetic_distance_mat_d2 = dist.nodes(tree_d2) ## Pairwise time of differences inferred with BEAST (ie taking into account substitution model)
colnames(genetic_distance_mat_d2) = c(names_seqs_d2, length(names_seqs_d2)+(1:(length(names_seqs_d2)-1)))
rownames(genetic_distance_mat_d2) = c(names_seqs_d2, length(names_seqs_d2)+(1:(length(names_seqs_d2)-1)))
hamming_dist_raw_d2 = genetic_distance_mat_d2

# D3
genetic_distance_mat_d3 = dist.nodes(tree_d3) ## Pairwise time of differences inferred with BEAST (ie taking into account substitution model)
colnames(genetic_distance_mat_d3) = c(names_seqs_d3, length(names_seqs_d3)+(1:(length(names_seqs_d3)-1)))
rownames(genetic_distance_mat_d3) = c(names_seqs_d3, length(names_seqs_d3)+(1:(length(names_seqs_d3)-1)))
hamming_dist_raw_d3 = genetic_distance_mat_d3

# D4
genetic_distance_mat_d4 = dist.nodes(tree_d4) ## Pairwise time of differences inferred with BEAST (ie taking into account substitution model)
colnames(genetic_distance_mat_d4) = c(names_seqs_d4, length(names_seqs_d4)+(1:(length(names_seqs_d4)-1)))
rownames(genetic_distance_mat_d4) = c(names_seqs_d4, length(names_seqs_d4)+(1:(length(names_seqs_d4)-1)))
hamming_dist_raw_d4 = genetic_distance_mat_d4

## Get the time each node
# D1
nroot_d1 = length(tree_d1$tip.label)+1 ## Root
distance_to_root_d1 = genetic_distance_mat_d1[nroot_d1,]
root_height_d1 = dataset_tips_d1$time[which(dataset_tips_d1$name_seq == names(distance_to_root_d1[4]))] - distance_to_root_d1[4]  ## Take one tip, doesn't matter which tip is used
nodes_height_d1 = root_height_d1 + distance_to_root_d1[length(names_seqs_d1)+(1:(length(names_seqs_d1)-1))]

# D2
nroot_d2 = length(tree_d2$tip.label)+1 ## Root
distance_to_root_d2 = genetic_distance_mat_d2[nroot_d2,]
root_height_d2 = dataset_tips_d2$time[which(dataset_tips_d2$name_seq == names(distance_to_root_d2[4]))] - distance_to_root_d2[4]  ## Take one tip, doesn't matter which tip is used
nodes_height_d2 = root_height_d2 + distance_to_root_d2[length(names_seqs_d2)+(1:(length(names_seqs_d2)-1))]

# D3
nroot_d3 = length(tree_d3$tip.label)+1 ## Root
distance_to_root_d3 = genetic_distance_mat_d3[nroot_d3,]
root_height_d3 = dataset_tips_d3$time[which(dataset_tips_d3$name_seq == names(distance_to_root_d3[4]))] - distance_to_root_d3[4]  ## Take one tip, doesn't matter which tip is used
nodes_height_d3 = root_height_d3 + distance_to_root_d3[length(names_seqs_d3)+(1:(length(names_seqs_d3)-1))]

# D4
nroot_d4 = length(tree_d4$tip.label)+1 ## Root
distance_to_root_d4 = genetic_distance_mat_d4[nroot_d4,]
root_height_d4 = dataset_tips_d4$time[which(dataset_tips_d4$name_seq == names(distance_to_root_d4[4]))] - distance_to_root_d4[4]  ## Take one tip, doesn't matter which tip is used
nodes_height_d4 = root_height_d4 + distance_to_root_d4[length(names_seqs_d4)+(1:(length(names_seqs_d4)-1))]
########################################################################################################################################

########################################################################################################################################
## Preparation data nodes
########################################################################################################################################
# Meta-data with nodes, for each serotype
dataset_with_nodes_d1 = data.frame('ID' = c(1:length(names_seqs_d1), length(names_seqs_d1)+(1:(length(names_seqs_d1)-1))),
                                   'name_seq' = c(names_seqs_d1, length(names_seqs_d1)+(1:(length(names_seqs_d1)-1))),
                                   'time' = c(dataset_tips_d1$time, nodes_height_d1),
                                   'is.node' = c(rep('no', length(names_seqs_d1)), rep('yes', (length(names_seqs_d1)-1))))
dataset_with_nodes_d1$genotype = meta_data_d1$genotype[match(dataset_with_nodes_d1$name_seq, meta_data_d1$shortname)]

dataset_with_nodes_d2 = data.frame('ID' = c(1:length(names_seqs_d2), length(names_seqs_d2)+(1:(length(names_seqs_d2)-1))),
                                   'name_seq' = c(names_seqs_d2, length(names_seqs_d2)+(1:(length(names_seqs_d2)-1))),
                                   'time' = c(dataset_tips_d2$time, nodes_height_d2),
                                   'is.node' = c(rep('no', length(names_seqs_d2)), rep('yes', (length(names_seqs_d2)-1))))
dataset_with_nodes_d2$genotype = meta_data_d2$genotype[match(dataset_with_nodes_d2$name_seq, meta_data_d2$shortname)]

dataset_with_nodes_d3 = data.frame('ID' = c(1:length(names_seqs_d3), length(names_seqs_d3)+(1:(length(names_seqs_d3)-1))),
                                   'name_seq' = c(names_seqs_d3, length(names_seqs_d3)+(1:(length(names_seqs_d3)-1))),
                                   'time' = c(dataset_tips_d3$time, nodes_height_d3),
                                   'is.node' = c(rep('no', length(names_seqs_d3)), rep('yes', (length(names_seqs_d3)-1))))
dataset_with_nodes_d3$genotype = meta_data_d3$genotype[match(dataset_with_nodes_d3$name_seq, meta_data_d3$shortname)]

dataset_with_nodes_d4 = data.frame('ID' = c(1:length(names_seqs_d4), length(names_seqs_d4)+(1:(length(names_seqs_d4)-1))),
                                'name_seq' = c(names_seqs_d4, length(names_seqs_d4)+(1:(length(names_seqs_d4)-1))),
                                'time' = c(dataset_tips_d4$time, nodes_height_d4),
                                'is.node' = c(rep('no', length(names_seqs_d4)), rep('yes', (length(names_seqs_d4)-1))))
dataset_with_nodes_d4$genotype = meta_data_d4$genotype[match(dataset_with_nodes_d4$name_seq, meta_data_d4$shortname)]
########################################################################################################################################

########################################################################################################################################
## Compute THD of every tip and node
########################################################################################################################################
## Window of time on which to search for samples in the population
wind = 1 ## years

bandwidth = thd.bandwidth(timescale, genome_length, mu_d1, q = 0.5) ## Corresponding bandwidth
dataset_with_nodes_d1$THD = THD_from_tips_and_nodes_time_windows(hamming_dist_raw_d1, tree_d1, wind, dataset_with_nodes_d1, mu = mu_d1)

bandwidth = thd.bandwidth(timescale, genome_length, mu_d2, q = 0.5) ## Corresponding bandwidth
dataset_with_nodes_d2$THD = THD_from_tips_and_nodes_time_windows(hamming_dist_raw_d2, tree_d2, wind, dataset_with_nodes_d2, mu = mu_d2)

bandwidth = thd.bandwidth(timescale, genome_length, mu_d3, q = 0.5) ## Corresponding bandwidth
dataset_with_nodes_d3$THD = THD_from_tips_and_nodes_time_windows(hamming_dist_raw_d3, tree_d3, wind, dataset_with_nodes_d3, mu = mu_d3)

bandwidth = thd.bandwidth(timescale, genome_length, mu_d4, q = 0.5) ## Corresponding bandwidth
dataset_with_nodes_d4$THD = THD_from_tips_and_nodes_time_windows(hamming_dist_raw_d4, tree_d4, wind, dataset_with_nodes_d4, mu = mu_d4)
########################################################################################################################################

########################################################################################################################################
## Generate colors needed
########################################################################################################################################
library(MetBrewer)
lev = levels(as.factor(c(dataset_with_nodes_d1$genotype, 
                          dataset_with_nodes_d2$genotype,
                          dataset_with_nodes_d3$genotype,
                          dataset_with_nodes_d4$genotype)))
n <- length(lev)
colors = met.brewer(name="Cross", n=n, type="continuous")

dataset_with_nodes_d1$genotype_color = dataset_with_nodes_d1$genotype
dataset_with_nodes_d1$genotype_color = as.factor(dataset_with_nodes_d1$genotype_color)
labels = levels(dataset_with_nodes_d1$genotype_color)
levels(dataset_with_nodes_d1$genotype_color) = colors[match(labels, lev)]
dataset_with_nodes_d1$genotype_color = as.character(dataset_with_nodes_d1$genotype_color)
dataset_with_nodes_d1$genotype_color[which(is.na(dataset_with_nodes_d1$genotype_color))] = 'grey60'

dataset_with_nodes_d2$genotype_color = dataset_with_nodes_d2$genotype
dataset_with_nodes_d2$genotype_color = as.factor(dataset_with_nodes_d2$genotype_color)
labels = levels(dataset_with_nodes_d2$genotype_color)
levels(dataset_with_nodes_d2$genotype_color) = colors[match(labels, lev)]
dataset_with_nodes_d2$genotype_color = as.character(dataset_with_nodes_d2$genotype_color)
dataset_with_nodes_d2$genotype_color[which(is.na(dataset_with_nodes_d2$genotype_color))] = 'grey60'

dataset_with_nodes_d3$genotype_color = dataset_with_nodes_d3$genotype
dataset_with_nodes_d3$genotype_color = as.factor(dataset_with_nodes_d3$genotype_color)
labels = levels(dataset_with_nodes_d3$genotype_color)
levels(dataset_with_nodes_d3$genotype_color) = colors[match(labels, lev)]
dataset_with_nodes_d3$genotype_color = as.character(dataset_with_nodes_d3$genotype_color)
dataset_with_nodes_d3$genotype_color[which(is.na(dataset_with_nodes_d3$genotype_color))] = 'grey60'

dataset_with_nodes_d4$genotype_color = dataset_with_nodes_d4$genotype
dataset_with_nodes_d4$genotype_color = as.factor(dataset_with_nodes_d4$genotype_color)
labels = levels(dataset_with_nodes_d4$genotype_color)
levels(dataset_with_nodes_d4$genotype_color) = colors[match(labels, lev)]
dataset_with_nodes_d4$genotype_color = as.character(dataset_with_nodes_d4$genotype_color)
dataset_with_nodes_d4$genotype_color[which(is.na(dataset_with_nodes_d4$genotype_color))] = 'grey60'
########################################################################################################################################

########################################################################################################################################
## Plot tree & THD below, with colors from clades
########################################################################################################################################
## Plot tree with reconstructed states
# pdf(file = 'Figure_raw_index_4DENVserotypes.pdf', width = 30/2.56, height = 8/2.56)
par(mfcol = c(2,4), oma = c(0,0,0,0), mar = c(4,4,0,0))
# dev.new()
min_year = 1940
max_year = 2020
max_index = 0.65

## Plot DENV 1
tree_d1 = ladderize(tree_d1, right = F)
plot(tree_d1, show.tip.label = FALSE, edge.width = 0.5, edge.color = 'grey', x.lim = c(min_year, max_year)-root_height_d1)
tiplabels(pch = 16, col = dataset_with_nodes_d1$genotype_color, cex = 0.5)
# nodelabels(pch = 16, col = dataset_with_nodes_d1$genotype_color[(length(dataset_with_nodes_d1$genotype_color)+1):length(dataset_with_nodes_d1$genotype_color)], cex = 0.5)
# axisPhylo(side = 1, root.time = root_height, backward = F)
axisPhylo_NL(side = 1, root.time = root_height_d1, backward = F,
             at_axis = seq(min_year, max_year, 10)-root_height_d1,
             lab_axis = seq(min_year, max_year, 10))
## Plot thd
plot(dataset_with_nodes_d1$time[which(dataset_with_nodes_d1$is.node == 'yes')], 
     dataset_with_nodes_d1$THD[which(dataset_with_nodes_d1$is.node == 'yes')], 
     col = adjustcolor(dataset_with_nodes_d1$genotype_color[which(dataset_with_nodes_d1$is.node == 'yes')], alpha.f = 1),
     bty = 'n', xlim = c(min_year, max_year), cex = 0.5,
     pch = 16, bty = 'n', ylim = c(0, max_index), 
     main = paste0(''), 
     ylab = 'Diversity index', xlab = 'Time (years)', yaxt = 'n', xaxt = 'n')
points(dataset_with_nodes_d1$time[which(dataset_with_nodes_d1$is.node == 'no')], 
       dataset_with_nodes_d1$THD[which(dataset_with_nodes_d1$is.node == 'no')], 
       col = adjustcolor(dataset_with_nodes_d1$genotype_color[which(dataset_with_nodes_d1$is.node == 'no')], alpha.f = 1),
       cex = 0.5, pch = 16)
axis(1, at = seq(min_year, max_year, 10), labels = seq(min_year, max_year, 10))
axis(2, las = 2)
legend('topright', 
       legend = dataset_with_nodes_d1$genotype[which(duplicated(dataset_with_nodes_d1$genotype_color)==F)],
       fill = dataset_with_nodes_d1$genotype_color[which(duplicated(dataset_with_nodes_d1$genotype_color)==F)],
       cex = 0.5, bty = 'n')

## Plot DENV 2
tree_d2 = ladderize(tree_d2, right = F)
plot(tree_d2, show.tip.label = FALSE, edge.width = 0.5, edge.color = 'grey', x.lim = c(min_year, max_year)-root_height_d2)
tiplabels(pch = 16, col = dataset_with_nodes_d2$genotype_color, cex = 0.5)
# nodelabels(pch = 16, col = dataset_with_nodes_d2$genotype_color[(length(dataset_with_nodes_d2$genotype_color)+1):length(dataset_with_nodes_d2$genotype_color)], cex = 1)
# axisPhylo(side = 1, root.time = root_height, backward = F)
axisPhylo_NL(side = 1, root.time = root_height_d2, backward = F,
             at_axis = seq(min_year, max_year, 10)-root_height_d2,
             lab_axis = seq(min_year, max_year, 10))
## Plot thd
plot(dataset_with_nodes_d2$time[which(dataset_with_nodes_d2$is.node == 'yes')], 
     dataset_with_nodes_d2$THD[which(dataset_with_nodes_d2$is.node == 'yes')], 
     col = adjustcolor(dataset_with_nodes_d2$genotype_color[which(dataset_with_nodes_d2$is.node == 'yes')], alpha.f = 1),
     bty = 'n', xlim = c(min_year, max_year), cex = 0.5,
     pch = 16, bty = 'n', ylim = c(0, max_index), 
     main = paste0(''), 
     ylab = 'Diversity index', xlab = 'Time (years)', yaxt = 'n', xaxt = 'n')
points(dataset_with_nodes_d2$time[which(dataset_with_nodes_d2$is.node == 'no')], 
       dataset_with_nodes_d2$THD[which(dataset_with_nodes_d2$is.node == 'no')], 
       col = adjustcolor(dataset_with_nodes_d2$genotype_color[which(dataset_with_nodes_d2$is.node == 'no')], alpha.f = 1),
       cex = 0.5, pch = 16)
axis(1, at = seq(min_year, max_year, 10), labels = seq(min_year, max_year, 10))
axis(2, las = 2)
legend('topright', 
       legend = dataset_with_nodes_d2$genotype[which(duplicated(dataset_with_nodes_d2$genotype_color)==F)],
       fill = dataset_with_nodes_d2$genotype_color[which(duplicated(dataset_with_nodes_d2$genotype_color)==F)],
       cex = 0.5, bty = 'n')

## Plot DENV 3
tree_d3 = ladderize(tree_d3, right = F)
plot(tree_d3, show.tip.label = FALSE, edge.width = 0.5, edge.color = 'grey', x.lim = c(min_year, max_year)-root_height_d3)
tiplabels(pch = 16, col = dataset_with_nodes_d3$genotype_color, cex = 0.5)
# nodelabels(pch = 16, col = dataset_with_nodes_d3$genotype_color[(length(dataset_with_nodes_d3$genotype_color)+1):length(dataset_with_nodes_d3$genotype_color)], cex = 1)
# axisPhylo(side = 1, root.time = root_height, backward = F)
axisPhylo_NL(side = 1, root.time = root_height_d3, backward = F,
             at_axis = seq(min_year, max_year, 10)-root_height_d3,
             lab_axis = seq(min_year, max_year, 10))
## Plot thd
plot(dataset_with_nodes_d3$time[which(dataset_with_nodes_d3$is.node == 'yes')], 
     dataset_with_nodes_d3$THD[which(dataset_with_nodes_d3$is.node == 'yes')], 
     col = adjustcolor(dataset_with_nodes_d3$genotype_color[which(dataset_with_nodes_d3$is.node == 'yes')], alpha.f = 1),
     bty = 'n', xlim = c(min_year, max_year), cex = 0.5,
     pch = 16, bty = 'n', ylim = c(0, max_index), 
     main = paste0(''), 
     ylab = 'Diversity index', xlab = 'Time (years)', yaxt = 'n', xaxt = 'n')
points(dataset_with_nodes_d3$time[which(dataset_with_nodes_d3$is.node == 'no')], 
       dataset_with_nodes_d3$THD[which(dataset_with_nodes_d3$is.node == 'no')], 
       col = adjustcolor(dataset_with_nodes_d3$genotype_color[which(dataset_with_nodes_d3$is.node == 'no')], alpha.f = 1),
       cex = 0.5, pch = 16)
axis(1, at = seq(min_year, max_year, 10), labels = seq(min_year, max_year, 10))
axis(2, las = 2)
legend('topright', 
       legend = dataset_with_nodes_d3$genotype[which(duplicated(dataset_with_nodes_d3$genotype_color)==F)],
       fill = dataset_with_nodes_d3$genotype_color[which(duplicated(dataset_with_nodes_d3$genotype_color)==F)],
       cex = 0.5, bty = 'n')

## Plot DENV 4
tree_d4 = ladderize(tree_d4, right = F)
plot(tree_d4, show.tip.label = FALSE, edge.width = 0.5, edge.color = 'grey', x.lim = c(min_year, max_year)-root_height_d4)
tiplabels(pch = 16, col = dataset_with_nodes_d4$genotype_color, cex = 0.5)
# nodelabels(pch = 16, col = dataset_with_nodes_d4$genotype_color[(length(dataset_with_nodes_d4$genotype_color)+1):length(dataset_with_nodes_d4$genotype_color)], cex = 1)
# axisPhylo(side = 1, root.time = root_height, backward = F)
axisPhylo_NL(side = 1, root.time = root_height_d4, backward = F,
             at_axis = seq(min_year, max_year, 10)-root_height_d4,
             lab_axis = seq(min_year, max_year, 10))
## Plot thd
plot(dataset_with_nodes_d4$time[which(dataset_with_nodes_d4$is.node == 'yes')], 
     dataset_with_nodes_d4$THD[which(dataset_with_nodes_d4$is.node == 'yes')], 
     col = adjustcolor(dataset_with_nodes_d4$genotype_color[which(dataset_with_nodes_d4$is.node == 'yes')], alpha.f = 1),
     bty = 'n', xlim = c(min_year, max_year), cex = 0.5,
     pch = 16, bty = 'n', ylim = c(0, max_index), 
     main = paste0(''), 
     ylab = 'Diversity index', xlab = 'Time (years)', yaxt = 'n', xaxt = 'n')
points(dataset_with_nodes_d4$time[which(dataset_with_nodes_d4$is.node == 'no')], 
       dataset_with_nodes_d4$THD[which(dataset_with_nodes_d4$is.node == 'no')], 
       col = adjustcolor(dataset_with_nodes_d4$genotype_color[which(dataset_with_nodes_d4$is.node == 'no')], alpha.f = 1),
       cex = 0.5, pch = 16)
axis(1, at = seq(min_year, max_year, 10), labels = seq(min_year, max_year, 10))
axis(2, las = 2)
legend('topright', 
       legend = dataset_with_nodes_d4$genotype[which(duplicated(dataset_with_nodes_d4$genotype_color)==F)],
       fill = dataset_with_nodes_d4$genotype_color[which(duplicated(dataset_with_nodes_d4$genotype_color)==F)],
       cex = 0.5, bty = 'n')
# dev.off()
########################################################################################################################################


