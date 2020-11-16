# Author: Samvida S. Venkatesh
# Date: 07/07/20

PATH = [redacted]

library(dplyr)

CHRS <- 1:22

BMI_score_files <- paste(PATH, "/BMI_GRS/BMI_GRS_CHR", 
                         CHRS, ".sscore", sep = "")

all_chrs <- lapply(CHRS, function (i) {
  data <- read.table(BMI_score_files[i], sep = "\t", header = T, 
                     comment.char = "?", stringsAsFactors = F)
  data$SCORE_SUM <- (data$NMISS_ALLELE_CT / 2)*(data$SCORE1_AVG)
  data$CHR <- i
  return (data[, c("IID", "NMISS_ALLELE_CT", "SCORE_SUM", "CHR")])
})
all_chrs <- bind_rows(all_chrs)

res <- all_chrs %>% group_by(IID) %>% summarise(nvars = sum(NMISS_ALLELE_CT / 2),
                                                FINAL_GRS = sum(SCORE_SUM) / nvars)

write.table(res[, c("IID", "FINAL_GRS")], 
            paste(PATH, "/results/BMI_GRS.txt", sep = ""),
            quote = F, row.names = F)