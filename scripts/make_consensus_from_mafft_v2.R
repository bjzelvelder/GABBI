#!/usr/bin/env Rscript
# Benjamin Zelvelder - 05/08/2025

# Check pairwise distance between sequences of an alignment and get a consensus sequence
# Usage: make_consensus_from_mafft.R $uce $rule $threshold 
# Example 1: Make majority rule consensus if sequences of uce-1002_pX_final.mafft.fasta have a pairwise distance of 95%
# make_consensus_from_mafft.R uce-1002_pX_final.mafft.fasta majority 0.05
# Example 2: Loop around all *.mafft.fasta sequences to make a iupac rule consensus if sequences of each fasta are 100% identical (only allowing length variation)
# make_consensus_from_mafft.R all iupac 0
# Example 3 : Majority rule consensus if less than 10% of sequences make up the second biggest cluster size (less stringent)
# make_consensus_from_mafft.R all majseq 0.05 0.1

library(ape)
library(stringr)
library(igraph)

# Def consensus functions
consensus_iupac_with_gaps <- function(alignment) {
  iupac_dict <- list(
    "a" = "a", "c" = "c", "g" = "g", "t" = "t",
    "ag" = "r", "ct" = "y", "cg" = "s", "at" = "w",
    "gt" = "k", "ac" = "m", "acg" = "v", "act" = "h",
    "agt" = "d", "cgt" = "b", "acgt" = "n",
    "-" = "-", "a-" = "a", "c-" = "c", "g-" = "g", "t-" = "t",
    "ag-" = "r", "ct-" = "y", "cg-" = "s", "at-" = "w",
    "gt-" = "k", "ac-" = "m", "acg-" = "v", "act-" = "h",
    "agt-" = "d", "cgt-" = "b", "acgt-" = "n"
  )
  # alignment = fasta
  seq_length <- ncol(alignment)
  consensus_seq <- character(seq_length)
  # i = 1
  for (i in 1:seq_length) {
    bases <- ape::base.freq(alignment[, i])
    bases_key <- paste0(sort(names(which(bases>0))), collapse = "")
    consensus_seq[i] <- ifelse(tolower(bases_key) %in% names(iupac_dict), iupac_dict[[bases_key]], "n")
  }
  return(consensus_seq)
}

consensus_majority <- function(alignment) {
  # alignment = fasta
  seq_length <- ncol(alignment)
  consensus_seq <- character(seq_length)
  # i = 605
  for (i in 1:seq_length) {
    bases <- ape::base.freq(alignment[, i])
    if (all(is.na(bases))) {
      consensus_seq[i] = "n"
    } else if (sum(bases==max(bases)) > 1) { # No base is majoritary, return "N"
      consensus_seq[i] = "n"
    } else {
      consensus_seq[i] = names(which(bases==max(bases)))
    }
  }
  return(consensus_seq)
}

args = commandArgs(trailingOnly=TRUE)
# setwd("~/Documents/these/Ochyromerini/Probe_design/IBA")
# args[1] = "all"
# args[1] = "uce-1002_pX_final.mafft.fasta"

cat("Working dir:",getwd(),"\n")

# First argument
uce = args[1]
if (uce == "all") {
  cat("Making consensus sequences of all alignments ending with .mafft.fasta.\n")
  files = list.files(pattern = "*.mafft.fasta$")
} else if (str_detect(uce,".fasta$") || str_detect(uce,".fas$")) {
  cat("Making consensus sequence of ",uce,".\n",sep="")
  files = uce
} else {
  stop("First argument must be either \"all\" or a specific .fasta file.", call.=FALSE)
}
# Second argument
# args[2] = "majseq"
cons_method = args[2]
if (cons_method == "iupac") {
  cat("Consensus with iupac method.\n")
} else if (cons_method == "majority" || cons_method == "majseq") {
  cat("Consensus with majority method.\n")
} else {
  stop("Second argument must be either \"iupac\", \"majseq\" or \"majority\".\n", call.=FALSE)
}

# Third argument
# dist_threshold = 0.05
dist_threshold = as.numeric(args[3])
if (!is.na(dist_threshold) && dist_threshold >= 0 && dist_threshold <= 1) {
  cat("Merging similar sequences if max pairwise genetic distance is smaller than ",dist_threshold,".\n")
} else {
  stop("Third argument must be a number between 0 and 1.\n", call.=FALSE)
}

# max_outlier_cluster_size = 0.2
if (cons_method == "majseq") {
  max_outlier_cluster_size = as.numeric(args[4])
  if (!is.na(max_outlier_cluster_size) && max_outlier_cluster_size >= 0 && max_outlier_cluster_size <= 1) {
    cat("Removing outliers if they make up less than",max_outlier_cluster_size*100,"% of sequences.\n")
  } else {
    stop("Fourth argument must be a number between 0 and 1.\n", call.=FALSE)
  }
}

# f = files[1]
# Loop over alignment(s) given
for (f in files) {
  # If all pairwise distance < dist_threshold, make consensus, otherwise, don't.
  fasta = read.dna(f,format = "fasta",as.matrix = T)
  outf = basename(gsub("_fasta","",gsub("\\.","_",f)))
  dist = as.matrix(dist.dna(fasta,model = "raw",pairwise.deletion = T))
  if (cons_method == "majseq") {
    # Find clusters of similar sequences with connected components
    dist = dist < dist_threshold
    dist[is.na(dist)] = FALSE # Change NA due to non overlapping sequences to FALSE (no reason to be connected)
    diag(dist) = FALSE
    g = graph_from_adjacency_matrix(dist, mode = "undirected", diag = FALSE)
    comps = components(g)
    # If the second biggest cluster is bigger than max_outlier_cluster_size, discard alignment, otherwise make majcons on the biggest cluster
    outlier_cluster = sort(comps$csize,decreasing = T)[2]
    if (!is.na(outlier_cluster) && outlier_cluster >= max_outlier_cluster_size*nrow(dist)) {
      cat(f,"discarded.",comps$no,"clusters of sequences detected exceeding the maximum outlier cluster size threshold.\n")
    } else {
      clusters = split(names(comps$membership), comps$membership)
      largest_cluster = clusters[[which.max(sapply(clusters, length))]]
      new_fasta = fasta[rownames(fasta) %in% largest_cluster,]
      consensus_seq = consensus_majority(new_fasta)
      consensus_seq_long = paste0(consensus_seq,collapse = "")
      header = paste0(">",outf," len=",length(consensus_seq))
      # Write consensus sequence
      out = file(paste0(f,".cons"))
      writeLines(c(header,consensus_seq_long),sep = "\n",out)
      close(out)
      cat(f,"consensus done.",if (!is.na(outlier_cluster)) "Outlier clusters removed.\n" else "\n")
    }
  } else {
    if (max(dist,na.rm = T) <= dist_threshold) {
      if (cons_method == "iupac") {
        consensus_seq = consensus_iupac_with_gaps(fasta)
      } else { # thus majority rule
        consensus_seq = consensus_majority(fasta)
      }
      consensus_seq_long = paste0(consensus_seq,collapse = "")
      header = paste0(">",outf," len=",length(consensus_seq))
      # Write consensus sequence
      out = file(paste0(f,".cons"))
      writeLines(c(header,consensus_seq_long),sep = "\n",out)
      close(out)
      cat(f,"consensus done.\n")
    } else {
      cat(f,"discarded. At least 1 sequence exceeded the maximum distance threshold.\n")
    }
  }
}

