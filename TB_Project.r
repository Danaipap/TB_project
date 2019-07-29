# Process of snp.tab files from each snippy output directory in R.
#  Set up working directory
setwd("path_to_snp.tab_files")

# read library readr
library(readr) 

# Due to large number of snp.tab files, we have processed initially each project accession number that we have retrieved from ENA separately. 
# READ & PROCESS SNP.TAB FILES FROM DATASET_1
# create a vector/list from all files
files <- list.files("path_to_snp.tab_files_of_each_project_accession_number") 

# Annotated with dataset extension
files <- paste("additional_path_to_dataset1_folder/", files, sep = "")  

# read data frames -Read snp.tab files from dataset1
allFiles.list <- lapply(files, read_delim, delim = '\t', col_names = TRUE)

# combine all data frames by combining the rows, remove duplicate headers
df <- do.call("rbind", allFiles.list) 

# create new dataframe with the columns you need- here we have kept column position, ref &alt ,Ftype, effect, locus_tag, gene and product
df = df[,c(2,4,5,7,11,12,13,14)] 

# remove vector/list after reading the initial dataframes
rm(allFiles.list,files) 

#  Process dataset1
#  load dplyr library
library(dplyr)

#  group  per position  & locustag
#  add an extra column, summarising the cumulative number of nucleotide variations in specific position
dataset1 <- df %>% group_by(POS, LOCUS_TAG) %>%  summarise(SNPs = n())

# remove interim dataframe
rm(df) 

#  Read & process the rest of datasets in the same way
#  We have read & processed the rest of the 16 datasets in the same way

# Further analysis of datasets(trim/bind/groupby)
# trim datasets#  for example for dataset1 and dataset2 see below
#  keep column position and summary of SNPs
data_a <- dataset1[,c(1,3)] 

#  keep column position and summary of SNPs
data_b <- dataset2[,c(1,3)] 

# Do the same for rest of the datasets

# Remove all non –trimmed datasets

rm(dataset1, dataset2, …) 

# rbind trimmed datasets
data <- rbind(data_a, data_b, …)

#  remove trimmed datasets
rm(data_a, data_b, …) 

# combined dataframes-group by POS and summarise nucleotide variants
library(dplyr) 
data <- data %>% group_by(POS) %>% 
  summarise(SNPs = sum(SNPs))

# load intersected file. The intersected file takes into account all the overlapping genomic positions in H37Rv and was created in bash shell with the command bedtools intersect

#  read & view intersected file 
intersected <- read_delim("~/TB_bioinformatics/data/intersected1","\t", escape_double = FALSE, trim_ws = TRUE)
view(intersected)

# remove the duplicate and empty columns from intersected file
intersected <- intersected[,c(1,3,5,6,7,8,9,10)]

# name columns in dataframe intersected
names(intersected) = c("GENOME","POS","START","STOP", "FTYPE", "LOCUS_TAG", "GENE","PRODUCT")

# create a dataframe with the same number of positions as H37Rv
wholegenome <- data.frame(POS = seq(1,4411532))

# merge intersected with wholegenome dataframe by the column position
intersected <- merge(intersected, wholegenome, by = "POS", all = TRUE)

# remove wholegenome dataframe# 
rm(wholegenome)

#  Remove NA values from all the columns of intersected dataframe (in this case it’s the non coding positions) & replace with appropriate names
intersected$GENOME[is.na(intersected$GENOME)] <- "NC_000962"
intersected$START[is.na(intersected$START)] <- "NON_CODING"
intersected$STOP[is.na(intersected$STOP)] <- "NON_CODING"
intersected$FTYPE[is.na(intersected$FTYPE)] <- "NON_CODING"
intersected$LOCUS_TAG[is.na(intersected$LOCUS_TAG)] <- "NON_CODING"
intersected$GENE[is.na(intersected$GENE)] <- "NON_CODING"
intersected$PRODUCT[is.na(intersected$PRODUCT)] <- "NON_CODING"

# Merge intersected dataframe with data  dataframe by column position
intersected_final <- merge(intersected, data, by = "POS", all = TRUE)

#  replace NA in SNP column with zero
intersected_final$SNPs[is.na(intersected_final$SNPs)] <- 0

# remove intersected dataframe
rm(intersected)

# create final dataframe, but exclude miscellaneous RNAs, mobile elements, repeat regions and non coding regions
df_NC <- intersected[!intersected$GENE=="NON_CODING",]
df_NO <- df_NC[!df_NC$FTYPE=="tRNA" 
&!df_NC$FTYPE=="repeat_region"&!df_NC$FTYPE=="mobile_element"&!df_NC$FTYPE=="ncRNA"&!df_NC$FTYPE=="misc_RNA"&!df_NC$FTYPE=="rRNA",]

# GROUP by Locus_tag to sum the snps per locus tag
df_final <- df_NO %>% group_by(LOCUS_TAG, GENE, PRODUCT,START, STOP) %>% #  group  per LOCUS_TAG # 
  summarise(SNPs = sum(SNPs))

# Normalise by gene length
# subset start –stop column
subfinal<- df_final[,c(4,5)] 

# make start and stop column numeric
subfinal$START <- as.numeric(as.character(subfinal$START)) 
subfinal$STOP <- as.numeric(as.character(subfinal$STOP))

# # subtract stop-start column to find every gene’s length
subfinal$LENGTH <- (subfinal$STOP-subfinal$START)

# merge subfinal with with df_final  dataframe
subfinal2 <- merge(df_final, subfinal, by = "START", all = TRUE)

# remove 2nd stop column
subfinal3 <- subfinal2[,c(1,2,3,4,5,6,8)]

# rename columns appropriately
names(subfinal3) <- c("START","LOCUS_TAG","GENE","PRODUCT","STOP","SNP","LENGTH")

# divide number of nucleotide variants/gene’s length to normalise value by gene’s length
final <- transform(subfinal3, RATIO= SNP / LENGTH) 

#  sort by lowest mutation/length ratio
final_sorted <- final[order(final$RATIO),]

# Summary statistics & 5th-95th quantile
summary(final_sorted)
quantile(final_sorted$RATE, probs = seq(0, 1, by= 0.05))

# Validation of results
# The files used for the validation were the files contained in each snippy directory containing the depth in each genomic position. 
# Each file was manipulated in bash shell and contained 1 column(depth column). 
# The first file in the directory used contained 2 columns (position and depth)

# set working directory
setwd("/path_to_files_containing_depth&position")

# create a list vector with file names
filenames <- list.files(full.names=TRUE)

# read all files from list
allFiles.list <- lapply(filenames, read.csv,header = TRUE)

# column bind all dataframes. 
# Note that the first dataframe on the list has 2 columns (position &depth).
# The rest have only one column(depth). The df dataframe will have the number of lines of H37Rv and 8536 columns(1st column with have the genomic positions)
df <- do.call("cbind", allFiles.list) 

# subset each gene on 5th percentile e.g csb gene from the large df file 
csb <- subset(df, POS > 967896 & POS < 968306)

# check the sum of genomic positions equal or above depth of 10- keep column 1(position) stable as first column
csb$above10 =apply(csb[,-1], 1, function(x) sum(abs(x) >= 10))

# keep datframe without the column with the sum >10 depth
csb_1 <- csb[,c(1:8536)]

# check the sum of genomic positions below depth of 10- keep column 1(position) stable as first column
csb_1$below10 = apply(csb_1[,-1], 1, function(x) sum(abs(x) < 10))

# rename csb to csb_new-it will contain >10 depth column
csb_new <-csb[,c(1,8537)]

# rename csb_1 dataframe to csb_1_new –will contain below 10 depth column
csb_1_new <- csb_1[,c(1,8537)]

# merge csb_new & csb_new_1 dataframes to create dataframe with >10 and <10 depth column
csb_merge <- merge(csb_1_new,csb_new, by="POS", all=TRUE)

# remove files you don’t needa
rm(csb_1,csb,csb_new , csb_1_new )

# save output
write.csv(csb_merge,"/reads/tsv_files/all_final/test/csb_merge", row.names = FALSE)
