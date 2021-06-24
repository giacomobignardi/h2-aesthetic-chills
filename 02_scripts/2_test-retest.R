#Author: Giacomo Bignardi
#Date: 21/03/2021
#
#
#
#Test-Retest Reliability Correlations
#Remove: remove pairs  with no survey completed
#Test-retest: calculate single item reliability score
#-Randomization: chose one member of the pair to avoid familiar resemblances' effects
#Program: Test-Retest-----------------------------------------------------------------------------------------------------------------------------

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
library(tidyverse)
library(tidylog)

twinWide = read.csv(sprintf("%s/%s/1_AesChillWideDZmf.csv", wdNOA,wdNOA_Data))[-1]
set.seed(42)
#_________________________________________________________________________________________________________
####Remove####
#remove participants without any information
twinWide = twinWide%>%filter(!(is.na(twinWide$neo43_7.1) & 
                                 is.na(twinWide$neo43_7.2) & 
                                 is.na(twinWide$neo43_8.1) & 
                                 is.na(twinWide$neo43_8.2) & 
                                 is.na(twinWide$neo43_10.1) &
                                 is.na(twinWide$neo43_10.2)))
#tidylog: filter: removed 320 rows (3%), 9,000 rows remaining

#save for later analysis
write_csv(twinWide,sprintf("%s/%s/2_AesChillWideDZmf.csv", wdNOA,wdNOA_Data))

####Test-retest####
####_Randomization----
####__Remove----
#drop twins for which less than one survey is avaiable
#index for complete scores 1:complete, 0:missings
twinWide$wave7_1  = ifelse(!is.na(twinWide$neo43_7.1),1,0)
twinWide$wave8_1  = ifelse(!is.na(twinWide$neo43_8.1),1,0)
twinWide$wave10_1 = ifelse(!is.na(twinWide$neo43_10.1),1,0)
twinWide$wave7_2  = ifelse(!is.na(twinWide$neo43_7.2),1,0)
twinWide$wave8_2  = ifelse(!is.na(twinWide$neo43_8.2),1,0)
twinWide$wave10_2 = ifelse(!is.na(twinWide$neo43_10.2),1,0)

#create index to filter twins with only one wave
twinWide$waves_1 = twinWide$wave7_1 + twinWide$wave8_1 + twinWide$wave10_1
twinWide$waves_2 =   twinWide$wave7_2 + twinWide$wave8_2 + twinWide$wave10_2
twinWide$waves = ifelse(twinWide$waves_1 <= 1 & twinWide$waves_2 <= 1 , 0, 1) #0: no more than one survey per twin pair is avaiable;
twinWide_testRetest = twinWide%>%filter(waves == 1)
#tidylog: removed 4,268 rows (47%), 4,732 rows remaining

#number of twins for which more than one survey is aviable
nrow(twinWide_testRetest%>%filter(waves_1>1)) + nrow(twinWide_testRetest%>%filter(waves_2>1))

#retain only b_df: pairs that have at least one twin with more than two surveys
####__First step: randomize for NA----
#take care of twins who have no complete answers, but that have a twin with at least one
#create index to point twins in one pair who have NA for all the survey
twinWide_testRetest$twinMissing = ifelse(is.na(twinWide_testRetest$neo43_7.1) & is.na(twinWide_testRetest$neo43_8.1) & is.na(twinWide_testRetest$neo43_10.1),2,
       ifelse(is.na(twinWide_testRetest$neo43_7.2) & is.na(twinWide_testRetest$neo43_8.2) & is.na(twinWide_testRetest$neo43_10.2),1,0))

#split dataframe for twins who have at least one complete score for both members of the pairs and who don't
twinWide_testRetest_both = twinWide_testRetest%>%filter(twinMissing ==0)
twinWide_testRetest_notBoth = twinWide_testRetest%>%filter(twinMissing !=0)

#create index to randomly select one of the to twin
twinWide_testRetest_both$twinWhich = floor(runif(nrow(twinWide_testRetest_both),min=1, max=3))
twinWide_testRetest_notBoth$twinWhich = twinWide_testRetest_notBoth$twinMissing
#bind back the two df
twinWide_testRetest = rbind(twinWide_testRetest_both,twinWide_testRetest_notBoth)


####__Second step: randomize for sex----
#Prepare dataframe for male
twinWide_testRetest_male = twinWide_testRetest%>%filter(twzyg == 1 | twzyg == 2 | twzyg == 5 | twzyg == 6)
#tidylog: removed 2,572 rows (54%), 2,160 rows remaining

#Prepare dataframe for female
twinWide_testRetest_female  = twinWide_testRetest%>%filter(twzyg == 3 | twzyg == 4 | twzyg == 5 | twzyg == 6)
#tidylog: removed 1,006 rows (21%), 3,726 rows remaining

#randomly select one of the two pairs based on the index
twinWide_testRetest_male_final = rbind(
  twinWide_testRetest_male[twinWide_testRetest_male$twinWhich ==1 & !twinWide_testRetest_male$twzyg == 6,
                         c("FID", "twzyg","neo43_7.1","neo43_8.1","neo43_10.1")]%>%
    rename(neo43_7 = "neo43_7.1",neo43_8 = "neo43_8.1",neo43_10 = "neo43_10.1"),
  twinWide_testRetest_male[twinWide_testRetest_male$twinWhich ==2 & !twinWide_testRetest_male$twzyg == 5,
                         c("FID", "twzyg","neo43_7.2","neo43_8.2","neo43_10.2")]%>%
    rename(neo43_7 = "neo43_7.2",neo43_8 = "neo43_8.2",neo43_10 = "neo43_10.2")
)
twinWide_testRetest_female_final = rbind(
  twinWide_testRetest_female[twinWide_testRetest_female$twinWhich ==1 & !twinWide_testRetest_female$twzyg == 5,
                            c("FID", "twzyg","neo43_7.1","neo43_8.1","neo43_10.1")]%>%
    rename(neo43_7 = "neo43_7.1",neo43_8 = "neo43_8.1",neo43_10 = "neo43_10.1"),
  twinWide_testRetest_female[twinWide_testRetest_female$twinWhich ==2 & !twinWide_testRetest_female$twzyg == 6,
                            c("FID", "twzyg","neo43_7.2","neo43_8.2","neo43_10.2")]%>%
    rename(neo43_7 = "neo43_7.2",neo43_8 = "neo43_8.2",neo43_10 = "neo43_10.2")
)

#sanity check (female + male == Total):
nrow(twinWide_testRetest_male_final) + nrow(twinWide_testRetest_female_final) == nrow(twinWide_testRetest)
twinWide_testRetest_final = rbind(twinWide_testRetest_male_final,twinWide_testRetest_female_final)

####Estimates----
testRetest_male = psych::corr.test(twinWide_testRetest_male_final[,3:5],use = "pairwise.complete.obs", adjust = "bonferroni")
testRetest_female = psych::corr.test(twinWide_testRetest_female_final[,3:5],use = "pairwise.complete.obs", adjust = "bonferroni")
print(testRetest_male, short = F)
print(testRetest_female, short = F)
