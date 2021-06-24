#Author: Giacomo Bignardi
#Date: 21/03/2021
#
#
#
#Randomized (optimized) imputation
#Aesthetic Chill = item 43 NEO-FFI Twin dataset arrange (Wide twin 1 and 2 in two column)
# After exclusion of twin pairs that did not report scores for any of the 3 surveys
# 1) we prioritized complete answers from twin pairs. For example, twin pairs were selected when both twins reported scores for the Item 43 on one of the surveys. 
#    If both twins answered in more than one survey, we randomly selected data for the pair from one of the complete surveys; 
#    + split 1
#    + split 2
#    + split 3
# 2) if one of the two twin’s response was missing for all of the three surveys, we randomly selected a survey with the response for the other twin
#    + split 1
#    + split 2
#    + split 3
#Program: AestheticChills_randomizedSelection-----------------------------------------------------------------------------------------------------------------------------
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

a_AesChill_DZmf_wide = read.csv(sprintf("%s/%s/2_AesChillWideDZmf.csv", wdNOA,wdNOA_Data))
#_________________________________________________________________________________________________________
#index for complete pair: 2 = both complete, 1 = one complete, 0 = both missings
a_AesChill_DZmf_wide$wave7  = ifelse(!is.na(a_AesChill_DZmf_wide$neo43_7.1) & !is.na(a_AesChill_DZmf_wide$neo43_7.2), 2,
                                   ifelse(!is.na(a_AesChill_DZmf_wide$neo43_7.1) | !is.na(a_AesChill_DZmf_wide$neo43_7.2), 1,
                                          0))
a_AesChill_DZmf_wide$wave8  = ifelse(!is.na(a_AesChill_DZmf_wide$neo43_8.1) & !is.na(a_AesChill_DZmf_wide$neo43_8.2), 2,
                                   ifelse(!is.na(a_AesChill_DZmf_wide$neo43_8.1) | !is.na(a_AesChill_DZmf_wide$neo43_8.2), 1,
                                          0))
a_AesChill_DZmf_wide$wave10 = ifelse(!is.na(a_AesChill_DZmf_wide$neo43_10.1) & !is.na(a_AesChill_DZmf_wide$neo43_10.2), 2,
                                   ifelse(!is.na(a_AesChill_DZmf_wide$neo43_10.1) | !is.na(a_AesChill_DZmf_wide$neo43_10.2), 1,
                                          0))

####COMP PAIRS####
####_first split: all pairs complete####
#Create matrix with all complete Complete pairs (ab)
a11_AesChill_DZmf_wide = a_AesChill_DZmf_wide%>%filter(!is.na(a_AesChill_DZmf_wide$neo43_7.1) & 
                                                            !is.na(a_AesChill_DZmf_wide$neo43_7.2) & 
                                                            !is.na(a_AesChill_DZmf_wide$neo43_8.1) & 
                                                            !is.na(a_AesChill_DZmf_wide$neo43_8.2) & 
                                                            !is.na(a_AesChill_DZmf_wide$neo43_10.1) &
                                                            !is.na(a_AesChill_DZmf_wide$neo43_10.2))
#tidylog: filter: removed 8,369 rows (93%), 631 rows remaining
#Everything else (ac)
a12_AesChill_DZmf_wide = a_AesChill_DZmf_wide%>%filter(!(!is.na(a_AesChill_DZmf_wide$neo43_7.1) & 
                                                             !is.na(a_AesChill_DZmf_wide$neo43_7.2) & 
                                                             !is.na(a_AesChill_DZmf_wide$neo43_8.1) & 
                                                             !is.na(a_AesChill_DZmf_wide$neo43_8.2) & 
                                                             !is.na(a_AesChill_DZmf_wide$neo43_10.1) &
                                                             !is.na(a_AesChill_DZmf_wide$neo43_10.2)))
#tidylog:filter: removed 631 rows (7%), 8,369 rows remaining



#random selection
set.seed(42)
#create a random index
a11_AesChill_DZmf_wide$rn = floor(runif(nrow(a11_AesChill_DZmf_wide), min=1, max=4))
#complete
a11_f_AesChill_DZmf_wide = rbind(
  a11_AesChill_DZmf_wide[a11_AesChill_DZmf_wide$rn ==1,c("FID", "twzyg","sex.1","age7.1","age7.2","neo43_7.1","neo43_7.2")]%>%
  rename(age.1 = "age7.1",age.2 = "age7.2",neo43.1 = "neo43_7.1",neo43.2 = "neo43_7.2"),
  a11_AesChill_DZmf_wide[a11_AesChill_DZmf_wide$rn ==2,c("FID", "twzyg","sex.1","age8.1","age8.2","neo43_8.1","neo43_8.2")]%>%
  rename(age.1 = "age8.1",age.2 = "age8.2",neo43.1 = "neo43_8.1",neo43.2 = "neo43_8.2"),
  a11_AesChill_DZmf_wide[a11_AesChill_DZmf_wide$rn ==3,c("FID", "twzyg","sex.1","age10.1","age10.2","neo43_10.1","neo43_10.2")]%>%
  rename(age.1 = "age10.1",age.2 = "age10.2",neo43.1 = "neo43_10.1",neo43.2 = "neo43_10.2")
)

####_second split: two pairs complete####
#Create matrix with  2 complete Complete pairs (acd)
a21_AesChill_DZmf_wide = a12_AesChill_DZmf_wide%>%filter((wave7==2 & wave8==2)|(wave7==2 & wave10==2)|(wave8==2 & wave10==2))
#tidylog: removed 6,897 rows (82%), 1,472 rows remaining
#Everything else (ace)
a22_AesChill_DZmf_wide = a12_AesChill_DZmf_wide%>%filter(!((wave7==2 & wave8==2)|(wave7==2 & wave10==2)|(wave8==2 & wave10==2)))
#tidylog: removed 1,472 rows (18%), 6,897 rows remaining

#random selection
set.seed(42)
#filter uncomplete pairs
a211_AesChill_DZmf_wide = a21_AesChill_DZmf_wide%>%filter(a21_AesChill_DZmf_wide$wave7 < 2)%>%select(-c("age7.1","age7.2","neo43_7.1","neo43_7.2","wave7"))
#create a random index
a211_AesChill_DZmf_wide$rn = floor(runif(nrow(a211_AesChill_DZmf_wide), min=1, max=3))

a212_AesChill_DZmf_wide = a21_AesChill_DZmf_wide%>%filter(a21_AesChill_DZmf_wide$wave8 < 2)%>%select(-c("age8.1","age8.2","neo43_8.1","neo43_8.2","wave8"))
#create a random index
a212_AesChill_DZmf_wide$rn = floor(runif(nrow(a212_AesChill_DZmf_wide), min=1, max=3))

a213_AesChill_DZmf_wide = a21_AesChill_DZmf_wide%>%filter(a21_AesChill_DZmf_wide$wave10 < 2)%>%select(-c("age10.1","age10.2","neo43_10.1","neo43_10.2","wave10"))
#create a random index
a213_AesChill_DZmf_wide$rn = floor(runif(nrow(a213_AesChill_DZmf_wide), min=1, max=3))


#complete
a21_f_AesChill_DZmf_wide = rbind(
  a211_AesChill_DZmf_wide[a211_AesChill_DZmf_wide$rn ==1,c("FID", "twzyg","sex.1","age8.1","age8.2","neo43_8.1","neo43_8.2")]%>%
    rename(age.1 = "age8.1",age.2 = "age8.2",neo43.1 = "neo43_8.1",neo43.2 = "neo43_8.2"),
  a211_AesChill_DZmf_wide[a211_AesChill_DZmf_wide$rn ==2,c("FID", "twzyg","sex.1","age10.1","age10.2","neo43_10.1","neo43_10.2")]%>%
    rename(age.1 = "age10.1",age.2 = "age10.2",neo43.1 = "neo43_10.1",neo43.2 = "neo43_10.2"),
  a212_AesChill_DZmf_wide[a212_AesChill_DZmf_wide$rn ==1,c("FID", "twzyg","sex.1","age7.1","age7.2","neo43_7.1","neo43_7.2")]%>%
    rename(age.1 = "age7.1",age.2 = "age7.2",neo43.1 = "neo43_7.1",neo43.2 = "neo43_7.2"),
  a212_AesChill_DZmf_wide[a212_AesChill_DZmf_wide$rn ==2,c("FID", "twzyg","sex.1","age10.1","age10.2","neo43_10.1","neo43_10.2")]%>%
    rename(age.1 = "age10.1",age.2 = "age10.2",neo43.1 = "neo43_10.1",neo43.2 = "neo43_10.2"),
  a213_AesChill_DZmf_wide[a213_AesChill_DZmf_wide$rn ==1,c("FID", "twzyg","sex.1","age7.1","age7.2","neo43_7.1","neo43_7.2")]%>%
    rename(age.1 = "age7.1",age.2 = "age7.2",neo43.1 = "neo43_7.1",neo43.2 = "neo43_7.2"),
  a213_AesChill_DZmf_wide[a213_AesChill_DZmf_wide$rn ==2,c("FID", "twzyg","sex.1","age8.1","age8.2","neo43_8.1","neo43_8.2")]%>%
    rename(age.1 = "age8.1",age.2 = "age8.2",neo43.1 = "neo43_8.1",neo43.2 = "neo43_8.2")
)

####_third split: one pair complete####
#Create matrix with one 1 complete Complete pairs (acd)
a31_AesChill_DZmf_wide = a22_AesChill_DZmf_wide%>%filter(wave7==2|wave8==2|wave10==2)
#tidylog: removed 3,865 rows (56%), 3,032 rows remaining
#Everything else (ace)
a32_AesChill_DZmf_wide = a22_AesChill_DZmf_wide%>%filter(!(wave7==2|wave8==2|wave10==2))
#tidylog: removed 3,032 rows (44%), 3,865 rows remaining

a31_f_AesChill_DZmf_wide = rbind(
  a31_AesChill_DZmf_wide[a31_AesChill_DZmf_wide$wave7 > 1,c("FID", "twzyg","sex.1","age7.1","age7.2","neo43_7.1","neo43_7.2")]%>%
  rename(age.1 = "age7.1",age.2 = "age7.2",neo43.1 = "neo43_7.1",neo43.2 = "neo43_7.2"),
  a31_AesChill_DZmf_wide[a31_AesChill_DZmf_wide$wave8 > 1,c("FID", "twzyg","sex.1","age8.1","age8.2","neo43_8.1","neo43_8.2")]%>%
  rename(age.1 = "age8.1",age.2 = "age8.2",neo43.1 = "neo43_8.1",neo43.2 = "neo43_8.2"),
  a31_AesChill_DZmf_wide[a31_AesChill_DZmf_wide$wave10 > 1,c("FID", "twzyg","sex.1","age10.1","age10.2","neo43_10.1","neo43_10.2")]%>%
  rename(age.1 = "age10.1",age.2 = "age10.2",neo43.1 = "neo43_10.1",neo43.2 = "neo43_10.2")
)

#End point 1
####UCOMP PAIRS####
#point 2

####_fourth split: 3 pairs uncomplete (at least one)####
#Create matrix with one 1 complete twin (acd)
a41_AesChill_DZmf_wide = a32_AesChill_DZmf_wide%>%filter(wave7== 1 & wave8==1 & wave10==1)
#filter: removed 3,562 rows (92%), 303 rows remaining
#Everything else (ace)
a42_AesChill_DZmf_wide = a32_AesChill_DZmf_wide%>%filter(!(wave7== 1 & wave8==1 & wave10==1))
#filter: removed 303 rows (8%), 3,562 rows remaining

#random selection
set.seed(42)
#create a random index
a41_AesChill_DZmf_wide$rn = floor(runif(nrow(a41_AesChill_DZmf_wide), min=1, max=4))
#complete
a41_f_AesChill_DZmf_wide = rbind(
  a41_AesChill_DZmf_wide[a41_AesChill_DZmf_wide$rn ==1,c("FID", "twzyg","sex.1","age7.1","age7.2","neo43_7.1","neo43_7.2")]%>%
    rename(age.1 = "age7.1",age.2 = "age7.2",neo43.1 = "neo43_7.1",neo43.2 = "neo43_7.2"),
  a41_AesChill_DZmf_wide[a41_AesChill_DZmf_wide$rn ==2,c("FID", "twzyg","sex.1","age8.1","age8.2","neo43_8.1","neo43_8.2")]%>%
    rename(age.1 = "age8.1",age.2 = "age8.2",neo43.1 = "neo43_8.1",neo43.2 = "neo43_8.2"),
  a41_AesChill_DZmf_wide[a41_AesChill_DZmf_wide$rn ==3,c("FID", "twzyg","sex.1","age10.1","age10.2","neo43_10.1","neo43_10.2")]%>%
    rename(age.1 = "age10.1",age.2 = "age10.2",neo43.1 = "neo43_10.1",neo43.2 = "neo43_10.2")
  )
  

####_fifth split: two pair incomplete (at least one) ########
a51_AesChill_DZmf_wide = a42_AesChill_DZmf_wide%>%filter((wave7==1 & wave8==1)|(wave7==1 & wave10==1)|(wave8==1 & wave10==1))
#tidylog:removed 2,514 rows (71%), 1,048 rows remaining
a52_AesChill_DZmf_wide = a42_AesChill_DZmf_wide%>%filter(!((wave7==1 & wave8==1)|(wave7==1 & wave10==1)|(wave8==1 & wave10==1)))
#tidylog:removed 1,048 rows (29%), 2,514 rows remaining

#random selection
set.seed(42)
#filter uncomplete pairs
a511_AesChill_DZmf_wide = a51_AesChill_DZmf_wide%>%filter(a51_AesChill_DZmf_wide$wave7 < 1)%>%select(-c("age7.1","age7.2","neo43_7.1","neo43_7.2","wave7"))
#create a random index
a511_AesChill_DZmf_wide$rn = floor(runif(nrow(a511_AesChill_DZmf_wide), min=1, max=3))

a512_AesChill_DZmf_wide = a51_AesChill_DZmf_wide%>%filter(a51_AesChill_DZmf_wide$wave8 < 1)%>%select(-c("age8.1","age8.2","neo43_8.1","neo43_8.2","wave8"))
#create a random index
a512_AesChill_DZmf_wide$rn = floor(runif(nrow(a512_AesChill_DZmf_wide), min=1, max=3))

a513_AesChill_DZmf_wide = a51_AesChill_DZmf_wide%>%filter(a51_AesChill_DZmf_wide$wave10 < 1)%>%select(-c("age10.1","age10.2","neo43_10.1","neo43_10.2","wave10"))
#create a random index
a513_AesChill_DZmf_wide$rn = floor(runif(nrow(a513_AesChill_DZmf_wide), min=1, max=3))


#complete
a51_f_AesChill_DZmf_wide = rbind(
  a511_AesChill_DZmf_wide[a511_AesChill_DZmf_wide$rn ==1,c("FID", "twzyg","sex.1","age8.1","age8.2","neo43_8.1","neo43_8.2")]%>%
    rename(age.1 = "age8.1",age.2 = "age8.2",neo43.1 = "neo43_8.1",neo43.2 = "neo43_8.2"),
  a511_AesChill_DZmf_wide[a511_AesChill_DZmf_wide$rn ==2,c("FID", "twzyg","sex.1","age10.1","age10.2","neo43_10.1","neo43_10.2")]%>%
    rename(age.1 = "age10.1",age.2 = "age10.2",neo43.1 = "neo43_10.1",neo43.2 = "neo43_10.2"),
  a512_AesChill_DZmf_wide[a512_AesChill_DZmf_wide$rn ==1,c("FID", "twzyg","sex.1","age7.1","age7.2","neo43_7.1","neo43_7.2")]%>%
    rename(age.1 = "age7.1",age.2 = "age7.2",neo43.1 = "neo43_7.1",neo43.2 = "neo43_7.2"),
  a512_AesChill_DZmf_wide[a512_AesChill_DZmf_wide$rn ==2,c("FID", "twzyg","sex.1","age10.1","age10.2","neo43_10.1","neo43_10.2")]%>%
    rename(age.1 = "age10.1",age.2 = "age10.2",neo43.1 = "neo43_10.1",neo43.2 = "neo43_10.2"),
  a513_AesChill_DZmf_wide[a513_AesChill_DZmf_wide$rn ==1,c("FID", "twzyg","sex.1","age7.1","age7.2","neo43_7.1","neo43_7.2")]%>%
    rename(age.1 = "age7.1",age.2 = "age7.2",neo43.1 = "neo43_7.1",neo43.2 = "neo43_7.2"),
  a513_AesChill_DZmf_wide[a513_AesChill_DZmf_wide$rn ==2,c("FID", "twzyg","sex.1","age8.1","age8.2","neo43_8.1","neo43_8.2")]%>%
    rename(age.1 = "age8.1",age.2 = "age8.2",neo43.1 = "neo43_8.1",neo43.2 = "neo43_8.2")
)



####_sixth split: only on pair incomplete (at least one)####
#Create matrix with one 1 complete Complete pairs (acd)
a61_AesChill_DZmf_wide = a52_AesChill_DZmf_wide%>%filter(wave7==1|wave8==1|wave10==1)#sanity check
#tidylog: no rows removed

a6_f_AesChill_DZmf_wide = rbind(
  a61_AesChill_DZmf_wide[a61_AesChill_DZmf_wide$wave7 > 0,c("FID", "twzyg","sex.1","age7.1","age7.2","neo43_7.1","neo43_7.2")]%>%
    rename(age.1 = "age7.1",age.2 = "age7.2",neo43.1 = "neo43_7.1",neo43.2 = "neo43_7.2"),
  a61_AesChill_DZmf_wide[a61_AesChill_DZmf_wide$wave8 > 0,c("FID", "twzyg","sex.1","age8.1","age8.2","neo43_8.1","neo43_8.2")]%>%
    rename(age.1 = "age8.1",age.2 = "age8.2",neo43.1 = "neo43_8.1",neo43.2 = "neo43_8.2"),
  a61_AesChill_DZmf_wide[a61_AesChill_DZmf_wide$wave10 > 0,c("FID", "twzyg","sex.1","age10.1","age10.2","neo43_10.1","neo43_10.2")]%>%
    rename(age.1 = "age10.1",age.2 = "age10.2",neo43.1 = "neo43_10.1",neo43.2 = "neo43_10.2")
)

####FINAL####
AesChillFinal <- rbind(a11_f_AesChill_DZmf_wide,a21_f_AesChill_DZmf_wide,a31_f_AesChill_DZmf_wide,a41_f_AesChill_DZmf_wide,a51_f_AesChill_DZmf_wide,a6_f_AesChill_DZmf_wide) #note that sex are missing for some of the uncomplete pairs

#additional exclusion: (age is missing)
a_AesChillFinal = AesChillFinal%>%filter(!(is.na(age.1) & is.na(age.2)))
#tidylog:removed 5 rows (<1%), 8,995 rows remaining

#deal with age
a11_AesChillFinal = a_AesChillFinal%>%filter(is.na(age.1) | is.na(age.2))
a11_AesChillFinal$rn = ifelse(is.na(a11_AesChillFinal$age.1),2,1)

a12_AesChillFinal = a_AesChillFinal%>%filter(!(is.na(age.1) | is.na(age.2)))
a12_AesChillFinal$rn = floor(runif(nrow(a12_AesChillFinal), min=1, max=3))

a11_12_AesChillFinal = rbind(a11_AesChillFinal,a12_AesChillFinal)
AesChillFinal_DOS = rbind(
  a11_12_AesChillFinal[a11_12_AesChillFinal$rn == 1,c("FID", "twzyg","age.1","neo43.1" ,"neo43.2")]%>%
    rename(age = "age.1", zyg = "twzyg", aesChill1 = "neo43.1" ,aesChill2 = "neo43.2"),
  a11_12_AesChillFinal[a11_12_AesChillFinal$rn == 2,c("FID", "twzyg","age.2","neo43.1" ,"neo43.2")]%>%
    rename(age = "age.2",zyg = "twzyg",aesChill1 = "neo43.1" ,aesChill2 = "neo43.2")
)


#Swap male female DOS
AesChillWideDos = rbind(AesChillFinal_DOS%>%filter(!zyg == 6),
                        AesChillFinal_DOS%>%filter(zyg == 6)%>%rename(aesChill1 = aesChill2, aesChill2 = aesChill1))

#collapse DOS male-female into cat 5
AesChillWideDos[which(AesChillWideDos$zyg == 6),"zyg"] = 5

#save_csv
write_csv(AesChillFinal_DOS,sprintf("%s/%s/3_AesChillWideDos_5-6.csv",wdNOA,wdNOA_Data))
write_csv(AesChillWideDos,sprintf("%s/%s/3_AesChillWideDos.csv",wdNOA,wdNOA_Data))


####DEMOGRAPHIC####
#age
AesChillWideDos%>%summarise(mean(age), sd(age), min(age), max(age))
#number of complete pairs
completePairs = AesChillWideDos%>%filter(!is.na(aesChill1) & !is.na(aesChill2))
#removed 3,863 rows (43%), 5,132 rows remaining
#frequencies of zygosity per twin pair
zygNum = as.data.frame(table(completePairs$zyg))
#MZ
zygNum[zygNum$Var1 == 1,  ]$Freq + zygNum[zygNum$Var1 == 3,  ]$Freq
#DZ
zygNum[zygNum$Var1 == 2,  ]$Freq + zygNum[zygNum$Var1 == 4,  ]$Freq
#DOS
zygNum[zygNum$Var1 == 5,  ]$Freq

#Combined sample [Table 1]
sexNum = as.data.frame(table(AesChillWideDos$zyg))
#MZ male
sexNum[sexNum$Var1 == 1,  ]$Freq
#MZ female
sexNum[sexNum$Var1 == 3,  ]$Freq 
#DZ male
sexNum[sexNum$Var1 == 2,  ]$Freq
#DZ female
sexNum[sexNum$Var1 == 4,  ]$Freq
#DOS
sexNum[sexNum$Var1 == 5,  ]$Freq
