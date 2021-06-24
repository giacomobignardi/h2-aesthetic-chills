# Program: oneACEvca.R
# Author: Giacomo Bignardi 
# adapted from:
# Hermine Maes (https://static1.squarespace.com/static/58b2481a9f7456906a3b9600/t/5e5b597055dd08365429864a/1583044977327/oneACEvca.pdf) and
# & Sarah Medland, Robert Kirkpatrick & Michael Hunter tutorial 2018 tutorial https://ibg.colorado.edu/cdrom2018/deleeuw/sexLimitation/sexLim2018.pdf
# Date: 21 03 2021
#
# Twin Univariate ACE model to estimate causes of variation across multiple groups
# Matrix style model - Raw data - Continuous data
####oneACEvca.R----------------------------------------------------------------------

#clear the Environment
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
library(OpenMx)
library(psych)
source(sprintf("%s/%s/functions/miFunctions.R",wdOA,wdOA_scripts))
# log:
# OpenMx version: 2.17.2 [GIT v2.17.2]
# R version: R version 3.6.2 (2019-12-12)
# Platform: x86_64-apple-darwin15.6.0 
# MacOS: 10.15.4
# Default optimizer: NPSOL
# NPSOL-enabled?: Yes
# OpenMP-enabled?: Yes
#load file
load(file = sprintf("%s/%s/4_SatModel_AesChill.RData",wdNOA,wdNOA_Data))

#PREPARE DATA ----------------------------------------------------------------------------------------------------------------------
#set seeds
set.seed(42)
# Load Data
colData <- dim(AesChillWideDos)
colData[2]
describe(AesChillWideDos[,1:colData[2]], skew=F)
# Select Variables for Analysis
vars <- 'AesChill' # list of variables names
nv <- 1 # number of variables
ntv <- nv*2 # number of total variables
selVars <- paste(vars,c(rep(1,nv),rep(2,nv)),sep="")

#rename data for the script
colnames(AesChillWideDos)[2] <- "zyg"
colnames(AesChillWideDos)[3] <- "age"
colnames(AesChillWideDos)[4] <- "AesChill1"
colnames(AesChillWideDos)[5] <- "AesChill2"

# Select Covariates for Analysis
covVars <- 'age'

# Select Data for Analysis
mzDatam <- subset(AesChillWideDos, zyg==1, c(selVars, covVars))
mzDataf <- subset(AesChillWideDos, zyg==3, c(selVars, covVars))
dzDatam <- subset(AesChillWideDos, zyg==2, c(selVars, covVars))
dzDataf <- subset(AesChillWideDos, zyg==4, c(selVars, covVars))
dosData <- subset(AesChillWideDos, zyg==5, c(selVars, covVars))

# Generate Descriptive Statistics
colMeans(mzDatam,na.rm=TRUE);
colMeans(mzDataf,na.rm=TRUE);
colMeans(dzDatam,na.rm=TRUE);
colMeans(dzDataf,na.rm=TRUE);
colMeans(dosData,na.rm=TRUE);
cov(mzDatam,use="complete")
cov(mzDataf,use="complete")
cov(dzDatam,use="complete")
cov(dzDataf,use="complete")
cov(dosData,use="complete")

# Set Starting Values
svBe <- 0.1 # start value for regressions
svMem <- 2.1 # start value for means
svMef <- 2.4 # start value for means
svPa <- .4 # start value for path coefficient
svPe <- .6 # start value for path coefficient 
#PREPARE MODEL----------------------------------------------------------------------------------------------------------------------

# Create Matrices for Covariates and linear Regression Coefficients
defL <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.age"), name="defL" )
pathBl <- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=svBe, label="b11", name="bl" )
# Create Algebra for expected Mean Matrices
meanGm <- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=svMem, labels="meanm", name="meanGm" )
meanGf <- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=svMef, labels="meanf", name="meanGf" )
# meanGo <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=c(svMem, svMef), labels=c("meanm","meanf"), name="meanGo")
expMeanMZm <- mxAlgebra( expression= cbind(meanGm + defL%*%bl,meanGm + defL%*%bl), name="expMeanMZm" )
expMeanDZm <- mxAlgebra( expression= cbind(meanGm + defL%*%bl,meanGm + defL%*%bl), name="expMeanDZm" )
expMeanMZf <- mxAlgebra( expression= cbind(meanGf + defL%*%bl,meanGf + defL%*%bl), name="expMeanMZf" )
expMeanDZf <- mxAlgebra( expression= cbind(meanGf + defL%*%bl,meanGf + defL%*%bl), name="expMeanDZf" )
expMeanDOS <- mxAlgebra( expression= cbind(meanGm + defL%*%bl,meanGf + defL%*%bl), name="expMeanDOS" )

# Create Matrices for Variance Components
covAm <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VA11m", name="VAm" )
covCm <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=0, label="VC11m", name="VCm" )
covEm <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPe, label="VE11m", name="VEm")
covAf <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VA11f", name="VAf" )
covCf <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=0, label="VC11f", name="VCf" )
covEf <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPe, label="VE11f", name="VEf" )

# Produce a vector which is = to 1 if variance is positive and -1 if variance is negative
signA <- mxAlgebra( ((-1)^omxLessThan(VAf,0))*((-1)^omxLessThan(VAm,0)), name="signA")
signC <- mxAlgebra( ((-1)^omxLessThan(VCf,0))*((-1)^omxLessThan(VCm,0)), name="signC")
# Create Matrix for Genetic Correlation between DOS twin pairs
rA  <- mxMatrix(type="Symm", nrow=nv, ncol=nv, free=TRUE, values=.5, label="rg",  name="rAos")

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covPm <- mxAlgebra( expression= VAm+VCm+VEm, name="Vm" )
covPf <- mxAlgebra( expression= VAf+VCf+VEf, name="Vf" )

# Define expected covariance
covMZm <- mxAlgebra( expression= VAm+VCm, name="cMZm" )
covDZm <- mxAlgebra( expression= 0.5%x%VAm+ VCm, name="cDZm" )
covMZf <- mxAlgebra( expression= VAf+VCf, name="cMZf" )
covDZf <- mxAlgebra( expression= 0.5%x%VAf+ VCf, name="cDZf" )
covDOS <- mxAlgebra( expression= rAos%x%((signA*sqrt(abs(VAf))%*%t(sqrt(abs(VAm))))+(signC*((sqrt(abs(VCf))%*%t(sqrt(abs(VCm))))))), name="cDOS")


#Create Algebra for expected Variance/Covariance Matrices in MZ, DZ and DOS twins
expCovMZm <- mxAlgebra( expression= rbind( cbind(Vm, cMZm), cbind(t(cMZm), Vm)), name="expCovMZm" )
expCovDZm <- mxAlgebra( expression= rbind( cbind(Vm, cDZm), cbind(t(cDZm), Vm)), name="expCovDZm" )
expCovMZf <- mxAlgebra( expression= rbind( cbind(Vf, cMZf), cbind(t(cMZf), Vf)), name="expCovMZf" )
expCovDZf <- mxAlgebra( expression= rbind( cbind(Vf, cDZf), cbind(t(cDZf), Vf)), name="expCovDZf" )
expCovDOS <- mxAlgebra( expression= rbind(cbind(Vm, cDOS),cbind(t(cDOS), Vf)),      name="expCovDOS" )

# Create Data Objects for Multiple Groups
dataMZm <- mxData( observed=mzDatam, type="raw" )
dataDZm <- mxData( observed=dzDatam, type="raw" )
dataMZf <- mxData( observed=mzDataf, type="raw" )
dataDZf <- mxData( observed=dzDataf, type="raw" )
dataDOS <- mxData( observed=dosData, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZm <- mxExpectationNormal( covariance="expCovMZm", means="expMeanMZm", dimnames=selVars )
expDZm <- mxExpectationNormal( covariance="expCovDZm", means="expMeanDZm", dimnames=selVars )
expMZf <- mxExpectationNormal( covariance="expCovMZf", means="expMeanMZf", dimnames=selVars )
expDZf <- mxExpectationNormal( covariance="expCovDZf", means="expMeanDZf", dimnames=selVars )
expDOS <- mxExpectationNormal( covariance="expCovDOS", means="expMeanDOS", dimnames=selVars )
funML <- mxFitFunctionML()

#Create Model Objects for Multiple Groups
#Create Model Objects for Multiple Groups
parsm <- list(pathBl, meanGm, covAm, covCm, covEm, covPm)
parsf <- list(pathBl, meanGf, covAf, covCf, covEf, covPf)
parsdos <- list(pathBl, meanGm, meanGf, covAm, covCm, covEm, covPm, covAf, covCf, covEf, covPf, covDOS, signA, signC, rA)
defs <- list( defL )

modelMZm <- mxModel( parsm, defs, expMeanMZm, covMZm, expCovMZm, dataMZm, expMZm, funML, name="MZm" )
modelDZm <- mxModel( parsm, defs, expMeanDZm, covDZm, expCovDZm, dataDZm, expDZm, funML, name="DZm" )
modelMZf <- mxModel( parsf, defs, expMeanMZf, covMZf, expCovMZf, dataMZf, expMZf, funML, name="MZf" )
modelDZf <- mxModel( parsf, defs, expMeanDZf, covDZf, expCovDZf, dataDZf, expDZf, funML, name="DZf" )
modelDOS <- mxModel( parsdos, defs, expMeanDOS, expCovDOS, dataDOS, expDOS, funML, name="DOS" )
multi <- mxFitFunctionMultigroup( c("MZm","DZm","MZf","DZf","DOS") )

#Create Algebra for Unstandardized and Standardized Variance Components
rowUSm <- rep('USm',nv)
rowUSf <- rep('USf',nv)
colUSm <- rep(c('VAm','VCm','VEm','SAm','SCm','SEm'),each=nv)
colUSf <- rep(c('VAf','VCf','VEf','SAf','SCf','SEf'),each=nv)
estUSm <- mxAlgebra( expression=cbind(VAm,VCm,VEm,VAm/Vm,VCm/Vm,VEm/Vm), name="USm", dimnames=list(rowUSm,colUSm) )
estUSf <- mxAlgebra( expression=cbind(VAf,VCf,VEf,VAf/Vf,VCf/Vf,VEf/Vf), name="USf", dimnames=list(rowUSf,colUSf) )

#Create Confidence Interval Objects
ciACEm <- mxCI( "USm[1,1:6]" )
ciACEf <- mxCI( "USf[1,1:6]" )

#Build Model with Confidence Intervals
modelACE <- mxModel( "oneACEvc", 
                     parsm, parsf, parsdos,
                     modelMZm, modelDZm,modelMZf, modelDZf, modelDOS, 
                     multi, 
                     estUSm, estUSf ,
                     ciACEf, ciACEm
                     )

#RUN MODEL----------------------------------------------------------------------------------------------------------------------

# Run ACE Model
mxCheckIdentification(modelACE, details=TRUE)
fitACE <- mxRun( modelACE, intervals= F)
sumACE <- summary( fitACE )
# Compare with Saturated Model

#TO RUN WITH SATURATED MODEL
#mxCompare( fitSAT, fitACE )

#if saturated model prior to genetic model
#lrtSAT(fitACE,4055.9346,1767)
# Print Goodness-of-fit Statistics & Parameter Estimates
fitGofs(fitACE)
fitEstCis(fitACE)
fitEsts(fitACE)

# RUN SUBMODELS----------------------------------------------------------------------------------------------------------------------
####_Sex limitiation Model----

#Expected Model AE
#Run AE model (male and female)
modelAE <- mxModel( fitACE, name="AE" )
modelAE <- omxSetParameters( modelAE, labels=c("VC11m","VC11f"), free=FALSE, values=0 )
fitAE <- mxRun( modelAE, intervals=F )
fitGofs(fitAE); fitEstCis(fitAE); fitEsts(fitAE)

#Model Comparision
#Run AE model: constrain mean equal for both sexes
modelAE_M<- mxModel( fitAE, name="AE_M" )
modelAE_M <- omxSetParameters( modelAE_M , label=c("meanm","meanf"), free=TRUE, 
                                   values= (svMem+svMem)/2  , newlabels='mean'  )
fitAE_M <- mxRun( modelAE_M , intervals=F )
fitGofs(fitAE_M); fitEstCis(fitAE_M); fitEsts(fitAE_M)

#Run Non-Scalar Sex-Limitation ACE model (= no qualitative sex differences) 
modelAE_Q <- mxModel( fitAE, name="AE_Q")
modelAE_Q <- omxSetParameters( modelAE_Q, labels="rg", free=FALSE, values=0.5)
fitAE_Q   <- mxRun( modelAE_Q)
fitGofs(fitAE_Q); fitEsts(fitAE_Q)

# Run AE model: constrain variance component to be equal for both sexes
modelAE_V <- mxModel( fitAE_Q, name="AE_V" )
modelAE_V <- omxSetParameters( modelAE_V, labels=c("VA11m","VA11f"), 
                                           free=TRUE, values= svPa, newlabels='VA11')
modelAE_V <- omxSetParameters( modelAE_V, labels=c("VE11m","VE11f"), 
                                           free=TRUE, values= svPe, newlabels='VE11')
fitAE_V <- mxRun( modelAE_V, intervals=T )
fitGofs(fitAE_V); fitEstCis(fitAE_V);fitEsts(fitAE_V)


####_compare Sex limitiation Model----
SexComp <- mxCompare( fitAE, subs <- list(fitAE_M, fitAE_Q, fitAE_V) )
mxCompare( fitSAT, fitAE)




####_Nested Variances SubModels----
####__Run ACE model----
#constrain variance component to be equal for both sexes
modelACE_V <- mxModel( fitACE, name="ACE_V" )
modelACE_V <- omxSetParameters( modelACE_V, labels="rg", free=FALSE, values=0.5)
modelACE_V <- omxSetParameters( modelACE_V, labels=c("VA11m","VA11f"), 
                                    free=TRUE, values= svPa, newlabels='VA11')
modelACE_V<- omxSetParameters( modelACE_V, labels=c("VC11m","VC11f"), 
                                    free=TRUE, values= 0, newlabels='VC11')
modelACE_V<- omxSetParameters( modelACE_V, labels=c("VE11m","VE11f"), 
                                    free=TRUE, values= svPe, newlabels='VE11')
fitACE_V <- mxRun( modelACE_V, intervals=F )
fitGofs(fitACE_V); fitEstCis(fitACE_V); fitEsts(fitACE_V)
####__Run CE model----
modelCE_V <- mxModel( fitACE_V, name="CE_V" )
modelCE_V <- omxSetParameters( modelCE_V, labels=c("VA11"), free=FALSE, values=0 )
fitCE_V <- mxRun( modelCE_V, intervals=F )
fitGofs(fitCE_V); fitEstCis(fitCE_V)
####__Run E model----
modelE_V <- mxModel( fitAE_V, name="E_V" )
modelE_V  <- omxSetParameters( modelE_V , labels=c("VA11"), free=FALSE, values=0 )
fitE_V <- mxRun( modelE_V , intervals=F )
fitGofs(fitE_V); fitEstCis(fitE_V)


####COMPARE MODEL###
# Print Comparative Fit Statistics
#mxCompare( fitSAT, fitACE_V)
ModCompSATandACE <- mxCompare(fitSAT, fitACE_V)
ModComp <- mxCompare( fitACE_V, nested <- list(fitAE_V, fitCE_V, fitE_V) )
VarianceNestedModel <- round(rbind(fitACE_V$USm$result,fitAE_V$USm$result,fitCE_V$USm$result,fitE_V$USm$result),4)

#estimates
fitACE_V$algebras$USm$result
fitAE_V$algebras$USm$result
fitCE_V$algebras$USm$result
fitE_V$algebras$USm$result
fitEstCis(fitAE_V)
fitEsts()

#SAVE OUTPUTS----------------------------------------------------------------------------------------------------------------------
save(VarianceNestedModel, SexComp, ModCompSATandACE, ModComp, modelAE_V, fitAE_V, file =sprintf("%s/%s/5_out_ACE.RData",wdOA,wdOA_output))

#sanity check: number of participants
mzDatam <- subset(mzDatam, zyg==1, c(selVars, covVars))
dzDatam <- subset(AesChillWideDos, zyg==2, c(selVars, covVars))
mzDataf <- subset(AesChillWideDos, zyg==3, c(selVars, covVars))
dzDataf <- subset(AesChillWideDos, zyg==4, c(selVars, covVars))
dosData <- subset(AesChillWideDos, zyg==5, c(selVars, covVars))
nrow(mzDatam)+nrow(dzDatam)+ nrow(dosData)/2 + nrow(mzDataf)+nrow(dzDataf)+ nrow(dosData)/2
