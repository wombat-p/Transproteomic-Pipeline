exp_design <-   read.csv("exp_design.txt", sep="\t")
exp_design <- exp_design[order(exp_design[,1]),,drop=F]
rownames(exp_design) <- exp_design[,1]
# adding 3rd columns with indices for replicate numbers
exp_design <- data.frame(exp_design, replicates=1)
counts <- vector("numeric", length(unique(exp_design[,2])))
names(counts) <- unique(exp_design[,2])
for (exp in 1:nrow(exp_design)) {
  counts[exp_design[exp, 2]] <- counts[exp_design[exp, 2]] + 1
  exp_design$replicates[exp] <- counts[exp_design[exp, 2]]
}  


## Protein table
# read and merge all output files from StPeter output converted to csv
quant_out <- quant_pep_out <-  list()
keep_columns_once <- c("protein_name","n_indistinguishable_proteins")
                       

keep_columns_all <- c("probability","total_number_peptides","total_number_distinct_peptides","Spectral.index",
                      "Normalized.spectral.index","number.of.quantified.peptides","Normalized.spectral.abundance.factor")
keep_pep_columns_once <- c("modified_peptide","peptide_sequence","charge","probability","calc_neutral_pep_mass",
                           "pseudo_name","subsuming_protein_entry")
keep_pep_columns_all <- c("fileName","protein_name","pct_spectrum_ids","confidence","initial_probability",
                          "nsp_adjusted_probability","fpkm_adjusted_probability","peptide_group_designator","weight",
                          "is_nondegenerate_evidence","n_instances","exp_tot_instances","is_contributing_evidence",
                          "StPeterQuant_peptide.sequence","StPeterQuant_peptide.charge","StPeterQuant_peptide.SI",
                          "StPeterQuant_peptide.SC")

# to map protein groups to their info
prot_info <- pep_info <- NULL
for (file in exp_design[,1]) {
  t_quant <- read.csv(sub(".raw", ".pep.interact.pep.prot_stpeter.prot_prot.csv", file))
  t_pep_quant <- read.csv(sub(".raw", ".pep.interact.pep.prot_stpeter.prot_pep.csv", file))
  # filter out non-quantified peptides
  t_pep_quant <- t_pep_quant[!is.na(t_pep_quant[,"StPeterQuant_peptide.sequence"]), ]
  
#  t_rownames <-  apply(t_quant[, c("protein_name","indistinguishable_proteins")], 1, paste, collapse=";")
  t_pep_rownames <-  apply(t_pep_quant[, c("modified_peptide","charge")], 1, paste, collapse=";")
  # some formatting of the secondary proteins as they appear in a json like format
#  t_rownames <- gsub("\\[|\\]|\\'", "", t_rownames)
#  rownames(t_quant) <- gsub(", ",";", t_rownames)
  rownames(t_pep_quant) <- t_pep_rownames
  rownames(t_quant) <- t_quant[,"protein_name"]
  t_prot_info <- t_quant[, keep_columns_once]
  t_pep_info <- t_pep_quant[, keep_pep_columns_once]
  quant_out[[file]] <- cbind(proteins=rownames(t_quant),t_quant[,keep_columns_all], stringsAsFactors=F)
  quant_pep_out[[file]] <- cbind(mod_peptides=rownames(t_pep_quant),t_pep_quant[,keep_pep_columns_all], stringsAsFactors=F)
  prot_info <- rbind(prot_info, t_prot_info)
  pep_info <- rbind(pep_info, t_pep_info)
  colnames(quant_out[[file]]) <- paste(colnames(quant_out[[file]]), exp_design[file,2], exp_design[file,3], sep="_")
  colnames(quant_pep_out[[file]]) <- paste(colnames(quant_pep_out[[file]]), exp_design[file,2], exp_design[file,3], sep="_")
}
prot_info <- unique(prot_info)
pep_info <- unique(pep_info)
all_quant <- Reduce(function(x, y) merge(x, y, by=1, all=TRUE), quant_out)
all_pep_quant <- Reduce(function(x, y) merge(x, y, by=1, all=TRUE), quant_pep_out)
all_quant <- cbind(prot_info[all_quant[,1],], all_quant[,2:ncol(all_quant)])
all_pep_quant <- cbind(pep_info[all_pep_quant[,1],], all_pep_quant[,2:ncol(all_pep_quant)])

# reorder columns
all_quant <- all_quant[,c(1,2,order(names(all_quant)[2:ncol(all_quant)]))]
all_pep_quant <- all_pep_quant[,c(1:7,order(names(all_pep_quant)[2:ncol(all_pep_quant)]))]

write.csv(all_quant, "all_prot_quant_merged.csv", row.names = F)
write.csv(all_pep_quant, "all_pep_quant_merged.csv", row.names = F)

## Peptide table (left out as tables contain PSMs and thus cannot be easily merged)
#quant_out <- list()
#keep_columns_once <- c("modified_peptide","protein","Quantification_SI","peptide_prev_aa","peptide_next_aa","modifications","peptide","num_missed_cleavages","calc_neutral_pep_mass")

#keep_columns_all <- c("spectrum","spectrumNativeID","assumed_charge","precursor_neutral_mass","retention_time_sec","start_scan","end_scan","xcorr","deltacn","deltacnstar","spscore","sprank","expect","num_matched_ions","tot_num_ions","massdiff","num_matched_peptides","fval","ntt","nmc","massd","peptideprophet_probability","peptideprophet_ntt_prob")

# to map protein groups to their info
#prot_info <- NULL
#for (file in exp_design[,1]) {
#  t_quant <- read.csv(sub(".raw", ".pep.interact.pep.prot_stpeter.prot_pep.csv", file, fixed = T))
#   rownames(t_quant) <-  apply(t_quant[, c("modified_peptide","assumed_charge")], 1 , paste, collapse=";")
#   print(keep_columns_once %in% colnames(t_quant))
#   print(keep_columns_all %in% colnames(t_quant))
#  t_prot_info <- t_quant[, keep_columns_once]
#  quant_out[[file]] <- t_quant[,keep_columns_all]
#  prot_info <- rbind(prot_info, t_prot_info)
#  colnames(quant_out[[file]]) <- paste(colnames(quant_out[[file]]), exp_design[file,2], exp_design[file,3], sep="_")
#}
#prot_info <- unique(prot_info)
#all_quant <- Reduce(function(x, y) merge(x, y, by=0, all=TRUE), quant_out)
#all_quant <- cbind(all_quant[,1], prot_info[all_quant[,1],], all_quant[,2:ncol(all_quant)])
#write.csv(all_quant, "all_pep_quant_merged.csv", row.names = F)

