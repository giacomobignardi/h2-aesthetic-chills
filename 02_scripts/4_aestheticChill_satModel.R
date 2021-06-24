# Program: oneSATca.R
# Author: Giacomo Bignardi adapted from Hermine Maes (https://static1.squarespace.com/static/58b2481a9f7456906a3b9600/t/5e5b57a80054c853617f6c79/1583044520354/oneSATca.pdf)
# Date: 21 03 2021
#
# Twin Univariate Saturated model to estimate means and (co)variances across multiple groups
# Matrix style model - Raw data - Continuous data
# oneSATca.R----------------------------------------------------------------------------------------------------------------------
#set seeds
set.seed(42)
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


library(OpenMx)
library(psych)
library(tidyverse)
source(sprintf("%s/functions/miFunctions.R",wdOA_scripts))

# log:
# OpenMx version: 2.17.2 [GIT v2.17.2]
# R version: R version 3.6.2 (2019-12-12)
# Platform: x86_64-apple-darwin15.6.0
# MacOS: 10.15.4
# Default optimizer: NPSOL
# NPSOL-enabled?: Yes
# OpenMP-enabled?: Yes

####PREPARE DATA----------------------------------------------------------------------------------
AesChillWideDos = read.csv(sprintf("%s/%s/3_AesChillWideDos.csv",wdNOA,wdNOA_Data))

colData = dim(AesChillWideDos)#
colData[2]
describe(AesChillWideDos[,1:colData[2]], skew=F)#
# Select Variables for Analysis
vars = 'aesChill' # list of variables names
nv = 1 # number of variables
ntv = nv*2 # number of total variables
selVars = paste(vars,c(rep(1,nv),rep(2,nv)),sep="")

#rename data for the script
colnames(AesChillWideDos)[2] = "zyg"
colnames(AesChillWideDos)[3] = "age"
colnames(AesChillWideDos)[4] = "aesChill1"
colnames(AesChillWideDos)[5] = "aesChill2"

# Select Covariates for Analysis
# AesChillWideDos[,'age'] = AesChillWideDos[,'age']/100
# AesChillWideDos = AesChillWideDos[-which(is.na(AesChillWideDos$age)),]
covVars = 'age'
# Select Data for Analysis
mzDatam = subset(AesChillWideDos, zyg==1, c(selVars, covVars))
mzDataf = subset(AesChillWideDos, zyg==3, c(selVars, covVars))
dzDatam = subset(AesChillWideDos, zyg==2, c(selVars, covVars))
dzDataf = subset(AesChillWideDos, zyg==4, c(selVars, covVars))
dosData = subset(AesChillWideDos, zyg==5, c(selVars, covVars))

# Generate Descriptive Statistics
colMeans(mzDatam,na.rm=TRUE);sd(mzDatam[,2])
colMeans(mzDataf,na.rm=TRUE);sd(mzDataf[,2])
colMeans(dzDatam,na.rm=TRUE);sd(dzDatam[,2])
colMeans(dzDataf,na.rm=TRUE);sd(dzDataf[,2])
colMeans(dosData,na.rm=TRUE);sd(dosData[,2])
cov(mzDatam,use="complete")
cov(mzDataf,use="complete")
cov(dzDatam,use="complete")
cov(dzDataf,use="complete")
cov(dosData,use="complete")
# Set Starting Values (play)
svBe = .1 # start value for regressions
svMem = 2.1 # start value for means male
svMef = 2.4 # start value for means female
svVa = 1.2 # start value for variance
lbVa = .1 # lower bound for variance

####PREPARE MODEL----------------------------------------------------------------------------------------------------------------------
# Create Matrices for Covariates and linear Regression Coefficients
defL = mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.age"), name="defL" )
pathBl = mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=svBe, label="b11", name="bl" )
# Create Algebra for expected Mean Matrices
meanMZm = mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMem, labels=c("mMZ1m","mMZ2m"), name="meanMZm" )
meanDZm = mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMem, labels=c("mDZ1m","mDZ2m"), name="meanDZm" )
meanMZf = mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMef, labels=c("mMZ1f","mMZ2f"), name="meanMZf" )
meanDZf = mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMef, labels=c("mDZ1f","mDZ2f"), name="meanDZf" )
meanDOS = mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=(svMem+svMef)/2, labels=c("mDOS1","mDOS2"), name="meanDOS" )
expMeanMZm = mxAlgebra( expression= meanMZm + cbind(defL%*%bl,defL%*%bl), name="expMeanMZm" )
expMeanDZm = mxAlgebra( expression= meanDZm + cbind(defL%*%bl,defL%*%bl), name="expMeanDZm" )
expMeanMZf = mxAlgebra( expression= meanMZf + cbind(defL%*%bl,defL%*%bl), name="expMeanMZf" )
expMeanDZf = mxAlgebra( expression= meanDZf + cbind(defL%*%bl,defL%*%bl), name="expMeanDZf" )
expMeanDOS = mxAlgebra( expression= meanDOS + cbind(defL%*%bl,defL%*%bl), name="expMeanDOS" )
# Create Algebra for expected Variance/Covariance Matrices
covMZm = mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv),
                    labels=c("vMZ1m","cMZ21m","vMZ2m"), name="covMZm" )
covDZm = mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv),
                    labels=c("vDZ1m","cDZ21m","vDZ2m"), name="covDZm" )
covMZf = mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv),
                    labels=c("vMZ1f","cMZ21f","vMZ2f"), name="covMZf" )
covDZf = mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv),
                    labels=c("vDZ1f","cDZ21f","vDZ2f"), name="covDZf" )
covDOS = mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv),
                    labels=c("vDOS1","cDOS21","vDOS2"), name="covDOS" )
# Calculate correlations from expected Covariance Matrices
corMZm = mxAlgebra( expression=cov2cor(covMZm), name ="corMZm")
corDZm = mxAlgebra( expression=cov2cor(covDZm), name ="corDZm")
corMZf = mxAlgebra( expression=cov2cor(covMZf), name ="corMZf")
corDZf = mxAlgebra( expression=cov2cor(covDZf), name ="corDZf")
corDOS =  mxAlgebra( expression=cov2cor(covDOS), name ="corDOS")
# Create Data Objects for Multiple Groups
dataMZm = mxData( observed=mzDatam, type="raw" )
dataDZm = mxData( observed=dzDatam, type="raw" )
dataMZf = mxData( observed=mzDataf, type="raw" )
dataDZf = mxData( observed=dzDataf, type="raw" )
dataDOS = mxData( observed=dosData, type="raw" )
# Create Expectation Objects for Multiple Groups
expMZm = mxExpectationNormal( covariance="covMZm", means="expMeanMZm", dimnames=selVars )
expDZm = mxExpectationNormal( covariance="covDZm", means="expMeanDZm", dimnames=selVars )
expMZf = mxExpectationNormal( covariance="covMZf", means="expMeanMZf", dimnames=selVars )
expDZf = mxExpectationNormal( covariance="covDZf", means="expMeanDZf", dimnames=selVars )
expDOS = mxExpectationNormal( covariance="covDOS", means="expMeanDOS", dimnames=selVars )
funML = mxFitFunctionML()
# Create Model Objects for Multiple Groups
pars = list( pathBl )
defs = list( defL )
modelMZm = mxModel( pars, defs, meanMZm, expMeanMZm, covMZm, corMZm, dataMZm, expMZm, funML, name="MZm" )
modelDZm = mxModel( pars, defs, meanDZm, expMeanDZm, covDZm, corDZm, dataDZm, expDZm, funML, name="DZm" )
modelMZf = mxModel( pars, defs, meanMZf, expMeanMZf, covMZf, corMZf, dataMZf, expMZf, funML, name="MZf" )
modelDZf = mxModel( pars, defs, meanDZf, expMeanDZf, covDZf, corDZf, dataDZf, expDZf, funML, name="DZf" )
modelDOS = mxModel( pars, defs, meanDOS, expMeanDOS, covDOS, corDOS, dataDOS, expDOS, funML, name="DOS" )
multi = mxFitFunctionMultigroup( c("MZm","DZm","MZf","DZf","DOS") )
# Create Confidence Interval Objects
ciCov = mxCI( c('MZm.covMZm','DZm.covDZm','MZf.covMZf','DZf.covDZf','DOS.covDOS') )
ciCor = mxCI( c('MZm.corMZm[2,1]','DZm.corDZm[2,1]','MZf.corMZf[2,1]','DZf.corDZf[2,1]','DOS.corDOS[2,1]') )
ciMean = mxCI( c('MZm.meanMZm','DZm.meanDZm', 'MZf.meanMZf','DZf.meanDZf', 'DOS.meanDOS') )
# Build Saturated Model with Confidence Intervals
modelSAT = mxModel( "oneSATca", pars, modelMZm, modelDZm, modelMZf, modelDZf, modelDOS, multi, ciCov,ciCor, ciMean)

####RUN MODEL----------------------------------------------------------------------------------------------------------------------
# Run Saturated Model
fitSAT = mxRun( modelSAT, intervals=F )
sumSAT = summary( fitSAT )
# Print Goodness-of-fit Statistics & Parameter Estimates
fitGofs( fitSAT )
fitEsts( fitSAT )
mxGetExpected( fitSAT, c("means","covariance") )

####RUN SUBMODEL----------------------------------------------------------------------------------------------------------------------
####_Covariate----
# Test Significance of Covariate
modelCOV = mxModel( fitSAT, name="COV" )
modelCOV = omxSetParameters( modelCOV, label="b11", free=FALSE, values=0 )
fitCOV = mxRun( modelCOV )

####_Birth Order(BO)----
# Constrain expected Means to be equal across Twin Order
modelBO_M = mxModel( fitSAT, name="BO_M" )
modelBO_M = omxSetParameters( modelBO_M, label=c("mMZ1m","mMZ2m"), free=TRUE, values=svMem, newlabels='mMZm' )
modelBO_M = omxSetParameters( modelBO_M, label=c("mDZ1m","mDZ2m"), free=TRUE, values=svMem, newlabels='mDZm' )
modelBO_M = omxSetParameters( modelBO_M, label=c("mMZ1f","mMZ2f"), free=TRUE, values=svMef, newlabels='mMZf' )
modelBO_M = omxSetParameters( modelBO_M, label=c("mDZ1f","mDZ2f"), free=TRUE, values=svMef, newlabels='mDZf' )
fitBO_M = mxRun( modelBO_M, intervals=F )
fitGofs(fitBO_M); fitEsts(fitBO_M)


#Constrain expected Means and Variances to be equal across Twin Order
modelBO_MV = mxModel( fitBO_M, name="BO_MV" )
modelBO_MV = omxSetParameters( modelBO_MV, label=c("vMZ1m","vMZ2m"), free=TRUE, values=svVa, newlabels='vMZm' )
modelBO_MV = omxSetParameters( modelBO_MV, label=c("vDZ1m","vDZ2m"), free=TRUE, values=svVa, newlabels='vDZm' )
modelBO_MV = omxSetParameters( modelBO_MV, label=c("vMZ1f","vMZ2f"), free=TRUE, values=svVa, newlabels='vMZf' )
modelBO_MV = omxSetParameters( modelBO_MV, label=c("vDZ1f","vDZ2f"), free=TRUE, values=svVa, newlabels='vDZf' )
fitBO_MV = mxRun( modelBO_MV, intervals=F)
fitGofs(fitBO_MV); fitEsts(fitBO_MV)


####_Zygosity----
# Constrain expected Means to be equal across Zygosity  (Taking Sex into account).
modelZyg_M = mxModel( fitBO_MV, name="Zyg_M" )
modelZyg_M  = omxSetParameters( modelZyg_M , label=c("mMZm","mDZm","mDOS1"), free=TRUE, values=(svMem+svMef)/2, newlabels='mm' )
modelZyg_M  = omxSetParameters( modelZyg_M , label=c("mMZf","mDZf","mDOS2"), free=TRUE, values=(svMem+svMef)/2, newlabels='mf' )
fitZyg_M = mxRun( modelZyg_M , intervals=F)
fitGofs(fitZyg_M); fitEsts(fitZyg_M)

# Constrain expected Mean and Variances to be equal across Zygosity (Taking Sex into account).
modelZyg_MV  = mxModel( fitZyg_M, name="Zyg_MV" )
modelZyg_MV  = omxSetParameters( modelZyg_MV , label=c("vMZm","vDZm","vDOS1"), free=TRUE, values= svVa, newlabels='vm' )
modelZyg_MV  = omxSetParameters( modelZyg_MV , label=c("vMZf","vDZf","vDOS2"), free=TRUE, values= svVa, newlabels='vf' )
fitZyg_MV = mxRun( modelZyg_MV , intervals=F)
fitGofs(fitZyg_MV); fitEsts(fitZyg_MV)


####_Sex----
# Constrain expected Means to be equal across Sexes.
modelSex_M = mxModel( fitZyg_MV, name="Sex_M" )
modelSex_M = omxSetParameters( modelSex_M, label=c('mm','mf'), free=TRUE, values=(svMem+svMef)/2, newlabels='m')
fitSex_M = mxRun( modelSex_M, intervals=F)
fitGofs(fitSex_M); fitEsts(fitSex_M)

# Constrain expected variances to be equal across Sexes.
modelSex_V = mxModel( fitZyg_MV, name="Sex_V" )
modelSex_V = omxSetParameters( modelSex_V, label=c('vm','vf'), free=TRUE, values= 1, newlabels='v')
fitSex_V = mxRun( modelSex_V, intervals=T)
fitGofs(fitSex_V); fitEsts(fitSex_V)

#Constrain expected covariances to be equal across MZ same sex.
modelSex_MZssCov = mxModel( fitSex_V, name="Sex_CMZ" )
modelSex_MZssCov = omxSetParameters( modelSex_MZssCov, label=c('cMZ21m','cMZ21f'), free=TRUE, values=.5, newlabels='cMZ')
fitSex_MZssCov = mxRun( modelSex_MZssCov, intervals=F)
fitGofs(fitSex_MZssCov); fitEsts(fitSex_MZssCov)


# Constrain expected covariances to be equal across DZ same sex.
modelSex_DZssCov = mxModel( fitSex_MZssCov, name="Sex_CDZss" )
modelSex_DZssCov = omxSetParameters( modelSex_DZssCov, label=c('cDZ21m','cDZ21f'), free=TRUE, values=.5, newlabels='cDZss')
fitSex_DZssCov = mxRun( modelSex_DZssCov, intervals=T)
fitGofs(fitSex_DZssCov); fitEsts(fitSex_DZssCov)

# Constrain expected covariances to be equal across DZ same sex and DOS (DZ).
modelSex_DOScov = mxModel( fitSex_DZssCov, name="Sex_CDZ" )
modelSex_DOScov = omxSetParameters( modelSex_DOScov, label=c('cDZss','cDOS21'), free=TRUE, values=.5, newlabels='cDZ')
fitSex_DOScov = mxRun(modelSex_DOScov, intervals=T)
fitGofs(fitSex_DOScov); fitEsts(fitSex_DOScov)


####COMPARE MODEL----
# Print Comparative Fit Statistics
SATComp = mxCompare(fitSAT, subs <- list(fitCOV, fitBO_M, fitBO_MV, fitZyg_M, fitZyg_MV, fitSex_M, fitSex_V, fitSex_MZssCov, fitSex_DZssCov, fitSex_DOScov) )

# SATcompTab = formattable::formattable(SATComp, align =c("l","c","c","c","c", "c", "c", "c", "r"))
####DEMOGRAPHICS####

#demographics extracted from the SAT
#Report sample size
#Check number of pairs cases 
fitGofS(fitSAT)
#log 
#Mx:oneSATca  
#statistics=14127  
#records=8995 
#parameters=26  
#constraints=0  
#df=14101  
#-2LL=43163.6023  
#cpu=4.3812
#optim=NPSOL  
#version=2.17.2  
#code=0
          
####_sample size----
#pairs
table(AesChillWideDos$zyg)
#complete pairs
AesChillWideDos$check = ifelse(!is.na(AesChillWideDos$aesChill1) & !is.na(AesChillWideDos$aesChill2),  1,  0 )
nrow(AesChillWideDos) - nrow(AesChillWideDos%>%filter(check == 0))
#complete pairs per zyg
#MZ
#complete pairs MZ male
nrow(AesChillWideDos[AesChillWideDos$zyg==1,]) - nrow(AesChillWideDos%>%filter(check == 0, zyg==1)) +
#complete pairs MZ female
nrow(AesChillWideDos[AesChillWideDos$zyg==3,]) - nrow(AesChillWideDos%>%filter(check == 0, zyg==3))
#DZ
#complete pairs DZ memale
nrow(AesChillWideDos[AesChillWideDos$zyg==2,]) - nrow(AesChillWideDos%>%filter(check == 0, zyg==2)) +
#complete pairs DZ female
nrow(AesChillWideDos[AesChillWideDos$zyg==4,]) - nrow(AesChillWideDos%>%filter(check == 0, zyg==4))
#DOS
#complete pairs oDZ 
nrow(AesChillWideDos[AesChillWideDos$zyg==5,]) - nrow(AesChillWideDos%>%filter(check == 0, zyg==5))

####_mean and standard deviation----
SATmodelAesChill = as.data.frame(fitEsts(fitSex_V))
sqrt(SATmodelAesChill[3,1])#sd
SATmodelAesChill[2,1]#mean male
SATmodelAesChill[6,1]#mean female

####_phenotypic correlation----
#Model 1: Fit Sex_V:the most parsimonious model in the hierarchical comparison within the saturated (sub)models with covarainces per sex
CisCOVmf = fitEstCis(fitSex_V)%>%as.data.frame()%>%rownames_to_column()

#Twin correlations
rmod1 = rbind(
MZmr = CisCOVmf[CisCOVmf$rowname == "MZm.corMZm[2,1]",2:4],
MZfr = CisCOVmf[CisCOVmf$rowname == "MZf.corMZf[2,1]",2:4],
DZmr = CisCOVmf[CisCOVmf$rowname == "DZm.corDZm[2,1]",2:4],
DZfr = CisCOVmf[CisCOVmf$rowname == "DZf.corDZf[2,1]",2:4],
DOSr = CisCOVmf[CisCOVmf$rowname ==  "DOS.corDOS[2,1]",2:4]
)%>%
  rownames_to_column()%>%
  rename(zygosity = "rowname")

#Model 2: fitSex_DZssCov:the most parsimonious model in the hierarchical comparison within the saturated (sub)models with covarainces per sex constrained within DZ same sex
CisCOVDOS = fitEstCis(fitSex_DZssCov)%>%as.data.frame()%>%rownames_to_column()

#Twin correlations
rmod2 = rbind(
  MZr = CisCOVDOS[CisCOVDOS$rowname == "MZm.corMZm[2,1]",2:4], #note that here m and f are equal
  DZssr = CisCOVDOS[CisCOVDOS$rowname == "DZm.corDZm[2,1]",2:4],
  DOSr = CisCOVDOS[CisCOVDOS$rowname == "DOS.corDOS[2,1]",2:4]
)%>%
  rownames_to_column()%>%
  rename(zygosity = "rowname")

#Model 3: fitSex_DOScov:the most parsimonious model in the hierarchical comparison within the saturated (sub)models with covarainces constrained within Zygositysame sex
CisCOV = fitEstCis(fitSex_DOScov)%>%as.data.frame()%>%rownames_to_column()

#Twin correlations
rmod3 = rbind(
  MZr = CisCOV[CisCOV$rowname == "MZm.corMZm[2,1]",2:4], #note that here m and f are equal
  DZr = CisCOV[CisCOV$rowname == "DZm.corDZm[2,1]",2:4] #note that here DZ and DOS are equal
)%>%
  rownames_to_column()%>%
  rename(zygosity = "rowname")

#save output enivronment 
save(fitSex_V, fitSex_DZssCov, fitSex_DOScov,fitSAT, AesChillWideDos, 
     SATComp,
     rmod1,rmod2,rmod3,
     file = sprintf("%s/%s/4_SatModel_AesChill.RData", wdNOA, wdNOA_Data))

#save ouput environemnt (OA)
save(fitSex_V, fitSex_DZssCov, fitSex_DOScov,fitSAT, 
     SATComp,
     rmod1,rmod2,rmod3,
     file = sprintf("%s/%s/4_out_SatModel_AesChill.RData", wdOA, wdOA_output))


