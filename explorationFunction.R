# BasicSummary function
BasicSummary <- function(df, dgts = 3){
  ## #
  ## ################################################################
  ## #
  ## # Create a basic summary of variables in the data frame df,
  ## # a data frame with one row for each column of df giving the
  ## # variable name, type, number of unique levels, the most
  ## # frequent level, its frequency and corresponding fraction of
  ## # records, the number of missing values and its corresponding
  ## # fraction of records
  ## #
  ## ################################################################
  ## #
  m <- ncol(df)
  varNames <- colnames(df)
  varType <- vector("character",m)
  topLevel <- vector("character",m)
  topCount <- vector("numeric",m)
  missCount <- vector("numeric",m)
  levels <- vector("numeric", m)
  
  for (i in 1:m){
    x <- df[,i]
    varType[i] <- class(x)
    xtab <- table(x, useNA = "ifany")
    levels[i] <- length(xtab)
    nums <- as.numeric(xtab)
    maxnum <- max(nums)
    topCount[i] <- maxnum
    maxIndex <- which.max(nums)
    lvls <- names(xtab)
    topLevel[i] <- lvls[maxIndex]
    missIndex <- which((is.na(x)) | (x == "") | (x == " "))
    missCount[i] <- length(missIndex)
  }
  n <- nrow(df)
  topFrac <- round(topCount/n, digits = dgts)
  missFrac <- round(missCount/n, digits = dgts)
  ## #
  summaryFrame <- data.frame(variable = varNames, type = varType,
                             levels = levels, topLevel = topLevel,
                             topCount = topCount, topFrac = topFrac,
                             missFreq = missCount, missFrac = missFrac)
  return(summaryFrame)
}

# Create a dataframe with numerical variables only
library(MASS)
df <- UScereal[c(2:8,10)]
# Defining the function as convert matrix 
convert_matrix <- function(df) {
  # I used sapply four times to calculate mean, variance, standard error and number of observations of the list elements using sapply. To covert to a matrix, I used rbind to bind the list row wise
  mytable <- rbind(sapply(df,mean,na.rm=TRUE),
                   sapply(df,median,na.rm=TRUE),
                   sapply(df,var,na.rm=TRUE),
                   sapply(df,sd,na.rm=TRUE),
                   sapply(df, function(df) max(df, na.rm=TRUE) - min(df, na.rm=TRUE)),
                   sapply(df, IQR, na.rm=TRUE))
  
  # Adding names to mytable
  dimnames(mytable)<-list(c("Mean","Median", "Variance","Standard Deviation","Range","Interquartile Range"),names(df))
  # Using round to get only two decimal places
  finaltable <- round(mytable,2)
  
  # returning mytable
  return(round(finaltable,2))  
}









#FindOutliers function

ThreeSigma <- function(x, t = 3){
  
  mu <- mean(x, na.rm = TRUE)
  sig <- sd(x, na.rm = TRUE)
  if (sig == 0){
    message("All non-missing x-values are identical")
  }
  up <- mu + t * sig
  down <- mu - t * sig
  out <- list(up = up, down = down)
  return(out)
}

Hampel <- function(x, t = 3){
  
  mu <- median(x, na.rm = TRUE)
  sig <- mad(x, na.rm = TRUE)
  if (sig == 0){
    message("Hampel identifer implosion: MAD scale estimate is zero")
  }
  up <- mu + t * sig
  down <- mu - t * sig
  out <- list(up = up, down = down)
  return(out)
}

BoxplotRule<- function(x, t = 1.5){
  
  xL <- quantile(x, na.rm = TRUE, probs = 0.25, names = FALSE)
  xU <- quantile(x, na.rm = TRUE, probs = 0.75, names = FALSE)
  Q <- xU - xL
  if (Q == 0){
    message("Boxplot rule implosion: interquartile distance is zero")
  }
  up <- xU + t * Q
  down <- xU - t * Q
  out <- list(up = up, down = down)
  return(out)
}   

ExtractDetails <- function(x, down, up){
  
  outClass <- rep("N", length(x))
  indexLo <- which(x < down)
  indexHi <- which(x > up)
  outClass[indexLo] <- "L"
  outClass[indexHi] <- "U"
  index <- union(indexLo, indexHi)
  values <- x[index]
  outClass <- outClass[index]
  nOut <- length(index)
  maxNom <- max(x[which(x <= up)])
  minNom <- min(x[which(x >= down)])
  outList <- list(nOut = nOut, lowLim = down,
                  upLim = up, minNom = minNom,
                  maxNom = maxNom, index = index,
                  values = values,
                  outClass = outClass)
  return(outList)
}
FindOutliers <- function(x, t3 = 3, tH = 3, tb = 1.5){
  threeLims <- ThreeSigma(x, t = t3)
  HampLims <- Hampel(x, t = tH)
  boxLims <- BoxplotRule(x, t = tb)
  
  n <- length(x)
  nMiss <- length(which(is.na(x)))
  
  threeList <- ExtractDetails(x, threeLims$down, threeLims$up)
  HampList <- ExtractDetails(x, HampLims$down, HampLims$up)
  boxList <- ExtractDetails(x, boxLims$down, boxLims$up)
  
  sumFrame <- data.frame(method = "ThreeSigma", n = n,
                         nMiss = nMiss, nOut = threeList$nOut,
                         lowLim = threeList$lowLim,
                         upLim = threeList$upLim,
                         minNom = threeList$minNom,
                         maxNom = threeList$maxNom)
  upFrame <- data.frame(method = "Hampel", n = n,
                        nMiss = nMiss, nOut = HampList$nOut,
                        lowLim = HampList$lowLim,
                        upLim = HampList$upLim,
                        minNom = HampList$minNom,
                        maxNom = HampList$maxNom)
  sumFrame <- rbind.data.frame(sumFrame, upFrame)
  upFrame <- data.frame(method = "BoxplotRule", n = n,
                        nMiss = nMiss, nOut = boxList$nOut,
                        lowLim = boxList$lowLim,
                        upLim = boxList$upLim,
                        minNom = boxList$minNom,
                        maxNom = boxList$maxNom)
  sumFrame <- rbind.data.frame(sumFrame, upFrame)
  
  threeFrame <- data.frame(index = threeList$index,
                           values = threeList$values,
                           type = threeList$outClass)
  HampFrame <- data.frame(index = HampList$index,
                          values = HampList$values,
                          type = HampList$outClass)
  boxFrame <- data.frame(index = boxList$index,
                         values = boxList$values,
                         type = boxList$outClass)
  outList <- list(summary = sumFrame, threeSigma = threeFrame,
                  Hampel = HampFrame, boxplotRule = boxFrame)
  return(outList)
}