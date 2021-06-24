#Author: Giacomo Bignardi
#Date: 28/01/2021
#
#
#
#tidy NTR df
#Aesthetic Chill = item 43 NEO-FFI Twin dataset arrange (Wide twin 1 and 2 in two column)
#Program: AestheticChills_tidy-----------------------------------------------------------------------------------------------------------------------------
#clear wd
rm(list = ls())

#set Open Access working directories
wdOA = getwd()
wdOA_scripts = "02_scripts"
wdOA_output = "03_outputs/processedData"

#set not Open Access working directories
wdNOA = substr(
  getwd(),
  0,
  nchar(getwd())-nchar("04_analysis_OA")-1
)
wdNOA_Data = "03_rawData/private"
wdNOA_ImageOutput = "05_images/image/processedData"

#load packages
library(haven) #required to read .sav 
library(dplyr)

#load twin data obtained from the NTR
dataChills = read_sav(sprintf("%s/%s/Giacomo\ Bignardi_Heritability\ of\ Aesthetic\ chills.sav", wdNOA,wdNOA_Data))
nrow(dataChills) #n cases 32781

#Select Twins 
AesChill = dataChills %>%
  dplyr::filter(Extension == 1 | Extension == 2) %>% #select twins
  dplyr::filter(multiple_type == 2) %>% #select doubles
  dplyr::select(c(FID, Extension, twzyg, sex, age7, age8, age10, neo43_7, neo43_8, neo43_10))%>%
  dplyr::filter(twzyg == 1 | twzyg == 2 | twzyg == 3 | twzyg == 4 | twzyg == 5 | twzyg == 6)


#Long to Wide
AesChillWideDZmf = reshape(as.data.frame(AesChill), timevar = "Extension" , idvar = c("FID","twzyg"), direction = "wide")
#save_csv
write.csv(AesChillWideDZmf,sprintf("%s/%s/1_AesChillWideDZmf.csv",wdNOA,wdNOA_Data))


