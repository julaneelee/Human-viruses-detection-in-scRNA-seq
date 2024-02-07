library(data.table)
library(dplyr)
library(stringr)

# Make the suite of efetch to fetch .fasta from NCBI
tbl <-  fread("path/to/Table_human_viruses.txt", header = FALSE)
tbl
tbl <- tbl[tbl$V6 != "Not available",]
tbl
nc_vec <-tbl[["V6"]]
nc_vec


# Loop over all possibilities
lapply(nc_vec, function(x){
  ncs <- str_trim(strsplit(x, ",")[[1]])
  lapply(ncs, function(one_nc){
    lapply(1:9, function(i){
      me <- paste0(one_nc, ".", i)
      data.frame(x = paste0("efetch -db nuccore -id ",me, ' -format fasta > ', me,'.fasta'))
    }) %>% rbindlist() 
  }) %>% rbindlist() 
}) %>% rbindlist() -> all_efetch
all_efetch

write.table(all_efetch, file = "step3_run_efetch.sh", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
