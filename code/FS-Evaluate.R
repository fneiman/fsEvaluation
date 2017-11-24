# fsEvaluation.R
# started: 11.23.2017
# R code based on Marko Procic 2013 "The goodness of fit and 
# and statstical significance of seriation solutions". Journal
# of Archaeological Science 40: 4552-4559.
# Fraser D. Neiman
# fneiman@monticello.org
# last upadte: 11.24.2017

library(ca)
library(plotrix)

setwd('C:\\Users\\fneiman\\Desktop\\fsEvaluation\\')

# read the data from the github repo
ezeroData <- read.csv('https://raw.github.com/fneiman/fsEvaluation/master/data/ezeroData.csv', sep= ',', stringsAsFactors = F)

# do the ca, sort the data on dim 1, and plot
dataMatrix <- as.matrix(ezeroData[,2:ncol(ezeroData)])
rownames(dataMatrix) <-  ezeroData[,1]
rowScores <- ca(dataMatrix)$rowcoord[,1]
sortedMatrix <- dataMatrix[order(rowScores),]
props <- prop.table(sortedMatrix,1) 
battleship.plot(props, col='grey')
# OK. now we have what we need: a data matrix and a set of scores to try out


# now we define three functions...

# define a function 
fsStats <- function(rowScores, dataMatrix){
# function to compute the number of Modes and the "Seriation Index"
# arguments:  rowscores - a vector of scores to order the rows (e.g. CA scores, MCDs)
#             dataMatrix - a matrix: counts of types(cols) in assemblages (rows).
# values:     nModes - the number of modes in the preder specified by the rowScores vector
#             sIndex - the Seriation Index (coefficient)
sortedMatrix <- dataMatrix[order(rowScores),]
props <- prop.table(sortedMatrix,1)
nRows <- nrow(props)
nCols <- ncol(props)
modeMat <- matrix(0, nRows,nCols) 
for (i in 2:(nRows-1)) {
  modeMat[i,] <- (props[i-1,] < props[i,]) & (props[i,] > props[i+1,])
}
modeMat[1,] <- (props[1,] > props[2,]) 
modeMat[nRows,] <- (props[nRows,] > props[nRows-1,]) 
nModes <- sum(modeMat)
maxModes <- ifelse(nRows %% 2 == 0, (nRows*nCols)/2, (nRows*(nCols+1))/2)
sIndex <- (maxModes -nModes)/ (maxModes - nCols)
list(nModes=nModes, sIndex=sIndex)
}

# Define a function 
fsStatsMC <- function(nPermutations, dataMatrix){
# function to compute the number of Modes and the "Seriation Index"
# for a given number of random permutations of a data matrix.  
# arguments:  nPermutations - number of permutations desired  
#             dataMatrix - a matrix: counts of types(cols) in assemblages (rows) 
  props <- prop.table(dataMatrix,1)
  nRows <- nrow(props)
  nCols <- ncol(props)
  mcValues <- matrix(0,nPermutations,2)
  colnames(mcValues) <- c('nModes', 'sIndex')
  for(k in 1:nPermutations){
    props <- props[order(runif(nRows)),]
    modeMat <- matrix(0, nRows,nCols)
    for (i in 2:(nRows-1)) {
      modeMat[i,] <- (props[i-1,] < props[i,]) & (props[i,] > props[i+1,])
    }
    modeMat[1,] <- (props[1,] > props[2,]) 
    modeMat[nRows,] <- (props[nRows,] > props[nRows-1,]) 
    nModes <- sum(modeMat)
    maxModes <- ifelse(nRows %% 2 == 0, (nRows*nCols)/2, (nRows*(nCols+1))/2)
    sIndex <- (maxModes -nModes)/ (maxModes - nCols)
    mcValues[k,] <- c(nModes, sIndex)
  }
mcValues
}


# define a function
fsAnalyze <- function(stats, statsMC ){
# function to compute MC probabilities, zscores, and z-score probabilities, based on
# the tesults from fsStats() and fsStatsMC()
# arguments:    StatsMC - the matrix of Monte Carlo values from fsStatsMC()
#               Stats  - the stats for the seriation solution from fsStats()
# empirical probabilty from the ECDF of the sIndex values
probMC_sIndex <- 1- (ecdf(statsMC[,2])(stats$sIndex))   
# z score for the Seriation Index. This is probably a good measure to compare 
# seriations on different datasets
z_sIndex <- (stats$sIndex - mean(statsMC[,2]))/ sd(statsMC[,2])
# prob(z >= z_sIndex)
prob_z_sIndex <- 1-pnorm(z_sIndex)
list('sIndex'= stats$sIndex,
      'probMC_sIndex' = probMC_sIndex, 
     'z_sIndex' = z_sIndex, 
     'prob_z_sIndex' = prob_z_sIndex )
}


# now we call the three functions and do a couple of plots
ezeroStats <- fsStats(rowScores, dataMatrix)
ezeroStats

ezeroStatsMC  <- fsStatsMC(1000, dataMatrix)

fsAnalyze (ezeroStats, ezeroStatsMC)

# first the Monte Carlo values for the number of modes to makes sure the result match 
# Porcic's figure 5. It does!

hist(ezeroStatsMC[,1], 20, col='grey',
     xlab= 'Number of Modes', ylab = 'Frequency')
abline(v=ezeroStats$nModes, col='blue', lwd=4)


# now do the analogous plots for the Seriaton Index values

hist(ezeroStatsMC[,2], 20, col='grey',
     xlab= 'Seriation Index', ylab = 'Frequency',
     xlim= c(0,1))
abline(v=ezeroStats$sIndex, col='blue', lwd=4)





