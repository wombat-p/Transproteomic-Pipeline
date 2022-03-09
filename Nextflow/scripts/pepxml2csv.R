######################################################################
# AUTH  Jeff Jones
# DATE  2018.07.23
# DESC  parse pepXML post TPP: xtandem > Tandem2XML > PeptideProphet
#
######################################################################

rm(list=ls())
require(XML)
library(progress)
options(warn=-1)

help_text <- "
 NAME
    pepxml2csv.R

 SYNOPSIS
    pepxml2csv.R --xml=<path_pepxml> --fdr=0.05 --rank=1

 DESCRIPTION
    extract peptide id data from PeptideProphet pepXML files for downstream analysis

 COMMAND LINE

    --xml <path_pepxml>

    --csv <path_csv> (optional: path_pepxml.csv)

    --fdr false discovery rate value between 0-1 (default: 0.01 ~ 1% FDR) 

    --rank max rank desired value >= 1 (default: 1) assumes user only wants highest ranking hit

 EXAMPLE

    Rscript pepxml2csv.R --xml=<path_to.pepXML>

"

###############################################################################
# USER INPUT
path_xml                      <- NULL
path_csv                      <- NULL
fdr_cutf                      <- 0.01
max_rank                      <- 1

for (arg in commandArgs()){
    arg_value <- as.character(sub("--[a-z]*\\=", "", arg))
    if( grepl("--xml", arg) ) path_xml <- arg_value
    if( grepl("--csv", arg) ) path_csv <- arg_value
    if( grepl("--fdr", arg) ) fdr_cutf <- arg_value
    if( grepl("--rank", arg) ) max_rank <- arg_value
    if( grepl("--help", arg) ) stop(help_text)
}

###############################################################################
# INPUT VALIDATION
message <- NULL
if(is.null(path_xml)) message <- stop("ERROR\n", "  no mzXML file declared\n")
if(!grepl("pep.*XML$", path_xml)) message <- paste0(message, "  mz file (--xml) not a supported format\n")
if(is.null(path_csv))
    path_csv = paste0(path_xml, ".csv")
if(!grepl(".csv$", path_csv)) message <- paste0(message, "  csv file (--csv) not a supported format\n")

if(!is.null(message)) stop("ERROR\n", message)

cat("pepXML to CSV started\n")
cat(" xml file:                        ", path_xml, "\n")
cat(" csv file:                        ", path_csv, "\n")
cat(" fdr cut off:                     ", fdr_cutf, "\n")
cat(" max rank:                        ", max_rank, "\n")

#
# Read in the data
#
cat(" reading xml file ...")
data <- xmlParse(path_xml)
file_xml <- sub('.*/', '', path_xml)

#
# convert to a master list
#
xml_data <- xmlToList(data)
cat("\n")


#
# probe for model performance values
#
roc_dat <- xml_data$analysis_summary$peptideprophet_summary$roc_error_data

if(is.null(roc_dat)) stop("ERROR\n", "  xml does not contain FDR statistics\n")

w_dp <- which(names(roc_dat) == 'roc_data_point')
w_ep <- which(names(roc_dat) == 'error_point')

roc_df <- do.call(rbind, lapply(lapply(roc_dat[w_dp], unlist), "[",
                                unique(unlist(c(sapply(roc_dat[w_dp],names))))))
roc_df <- as.data.frame(roc_df)
roc_df$analysts <- row.names(roc_df)

err_df <- do.call(rbind, lapply(lapply(roc_dat[w_ep], unlist), "[",
                                unique(unlist(c(sapply(roc_dat[w_ep],names))))))
err_df <- as.data.frame(err_df)
err_df$analysts <- row.names(err_df)
err_df$sensitivity <- NA

pdf <- rbind(roc_df, err_df)
row.names(pdf) <- 1:dim(pdf)[1]

pdf$min_prob <- as.numeric(as.character(pdf$min_prob))
pdf$sensitivity <- as.numeric(as.character(pdf$sensitivity))
pdf$error <- as.numeric(as.character(pdf$error))
pdf$num_corr <- as.numeric(as.character(pdf$num_corr))
pdf$num_incorr <- as.numeric(as.character(pdf$num_incorr))

w_pcf <- which(pdf$error == max(pdf[pdf$error <= fdr_cutf &
                                        pdf$analysts == 'error_point',]$error) &
                   pdf$analysts == 'error_point')
min_prob <- pdf[w_pcf,]$min_prob
cat("\n")
print(pdf[w_pcf,], row.names = FALSE)
cat("\n")


#
# traverse the master list to pull id data
#
w_specs <- which(names(xml_data$msms_run_summary) == 'spectrum_query')
pb <- progress_bar$new(
    format = " extracting data [:bar] :percent eta: :eta",
    total = length(w_specs), clear = FALSE, width= 60)
hits <- list()
for ( spec_i in w_specs ){
    
    pb$tick()
    
    hits[[spec_i]] <- list()
    
    spec <- xml_data$msms_run_summary[spec_i]
    
    #
    # grab the spectrum attributes
    #
    spec_att <- append(unlist(spec$spectrum_query$.attrs), file_xml)
    names(spec_att)[length(spec_att)] <- 'fileName'
    
    for(result in spec$spectrum_query$search_result){
        
        #
        # grab the result attributes
        # 
        result_att <- unlist(result$.attrs)
        
        #
        # grab the X!T scoring values 
        #    
        w_search_scores <- which(names(result) == 'search_score')
        score_name <- as.character(as.data.frame(t(as.data.frame(result[w_search_scores])))$name)
        score_vals <- as.character(as.data.frame(t(as.data.frame(result[w_search_scores])))$value)
        names(score_vals) <- score_name  
        
        #
        # grab the Prophet attributes
        #
        prophet_att <- unlist(result$analysis_result$peptideprophet_result$.attrs)
        
        #
        # grab the Prophet scores
        #    
        prophet_scores <- result$analysis_result$peptideprophet_result$search_score_summary
        w_prophet_param <- which(names(result) == 'parameter')
        prophet_name <- as.character(as.data.frame(t(as.data.frame(prophet_scores[w_prophet_param])))$name)
        prophet_vals <- as.character(as.data.frame(t(as.data.frame(prophet_scores[w_prophet_param])))$value)
        names(prophet_vals) <- prophet_name  
        
        #
        # append all the values to our new list
        #
        hits[[spec_i]] <- append(hits[[spec_i]], as.list(spec_att))
        hits[[spec_i]] <- append(hits[[spec_i]], as.list(result_att))
        hits[[spec_i]] <- append(hits[[spec_i]], as.list(score_vals))
        hits[[spec_i]] <- append(hits[[spec_i]], as.list(prophet_att))
        hits[[spec_i]] <- append(hits[[spec_i]], as.list(prophet_vals))
    }
}

#
# use lapply to converge the hits list into a matrix with named columns
#
data <- do.call(rbind, lapply(lapply(hits, unlist), "[",
                              unique(unlist(c(sapply(hits,names))))))

#
# convert the matrix to a data frame
#
data <- as.data.frame(data)
colnames(data) <- unique(unlist(c(sapply(hits,names))))
data$probability <- as.numeric(as.character(data$probability))
data$hit_rank <- as.numeric(as.character(data$hit_rank))


#
# filter data by cutoff min_prob and rank
#
data <- data[data$probability >= min_prob & data$hit_rank <= max_rank & !is.na(data$spectrum),]

write.csv(data, path_csv, row.names = FALSE)
cat("\nFINISHED:\n")
cat(" data written to                  ", path_csv, "\n")

path_csv <- sub(".csv", "_performance.csv", path_csv)
write.csv(pdf, path_csv, row.names = FALSE)
cat(" performance roc written to       ", path_csv, "\n")