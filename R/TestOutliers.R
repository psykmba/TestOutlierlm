if(!require(olsrr)){
  install.packages("olsrr")
  library(olsrr)
}
if(!require(EnvStats)){
  install.packages("EnvStats")
  library(EnvStats)
}

#' testOutliers
#' This is the main function for the package, it takes a lm formula
#' and data and handels outliers in different ways. Thereafter it
#' recalculates the lm model with data that has been handled
#' @param form an lm formula
#' @param data data that should include all variables in the lm formula
#' @param limit sets alpha that is used for some outlier detection models
#' @param missing not used
#' @param otherVar not used
#'
#' @return a testOutliers object that include summary of all models
#'        and all datafiles, including the original.
#' @export
testOutliers <- function(form, data, limit = .001, missing = "None", otherVar = c())
{
  getInsideRange <- function(s, r)
  {
    return (ifelse(s >= r[1] & s <= r[2], s, ifelse(s < r[1], r[1], r[2])))
  }
  deleteOutsideRange <- function( s, r)
  {
    return(ifelse (s < r[1], NA, ifelse(s > r[2], NA, s) ))
  }

  Mahalanobis <- function()
  {
    scaleCor <- stats::cov(data[,-1])
    Outliers <- stats::mahalanobis(data[,-1], colMeans(data[,-1]), scaleCor)
    out <-  Outliers > stats::qchisq(1-limit, length(data))
    newData <- data
    newData[out,] <- NA
    return(newData)
  }
  Leverage <- function(model)
  {
    Outliers <- stats::hatvalues(model)
    newData <- dplyr::filter(data, Outliers < (2*(ncol(data)))/nrow(data))
    out <-  Outliers > (2*(ncol(data)))/nrow(data)
    newData <- data
    newData[out,] <- NA
    return(newData)
  }

  SD <- function()
  {
    newFrame <-  data.frame(row.names = 1:nrow(data))
     for(scale in data)
    {
      m <- mean(scale)
      sd <- sd(scale) * stats::qnorm(1 - limit)
      r <- range(m+sd, m-sd)
       newFrame <- cbind(newFrame, deleteOutsideRange(scale, r))
    }
     names(newFrame) <- names(data)


    return(newFrame)
  }
  Change <- function()
  {
    newFrame <- data.frame(row.names = 1:nrow(data))
     for(scale in data)
    {
      m <- mean(scale)
      sd <- sd(scale) * stats::qnorm(1 - limit)
      r <- range(m+sd, m-sd)
      newFrame <- cbind(newFrame, getInsideRange(scale, r))
    }

    names(newFrame) <- names(data)
    View(newFrame)

    return(newFrame)
  }
  CooksDistance4N <- function(model)
  {
    cooksd <- cooks.distance(model)
    sample_size <- nrow(data)
    plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
    abline(h = 4/sample_size, col="red")  # add cutoff line
    text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4/sample_size, names(cooksd),""), col="red")  # add labels
    influential <- as.numeric(names(cooksd)[(cooksd > (4/sample_size))])
    newData <- data
    newData[influential, ] <- NA
     return(newData)

  }
  DFBetas <- function(model)
  {
    dfb <- dfbetas(model)
    n <- length(model$residuals)

    #calculate DFBETAS threshold value
    thresh <- 2/sqrt(n)
     del <- apply(dfb[,-1], MARGIN = 1,
                 FUN = function(x) {return(any(abs(x)>thresh))})
     newData <- data
     newData[!del, ] <- NA
     return(newData)

  }
  DFBetasMissing <- function(model)
  {
    dfb <- dfbetas(model)
    n <- length(model$residuals)

    #calculate DFBETAS threshold value
    thresh <- 2/sqrt(n)
    del <- apply(dfb[,-1], MARGIN = 2,
                 FUN = function(x) {return(ifelse(abs(x)>thresh, FALSE, TRUE))})
    tillf <- data[,-1]
     for (i in 1:ncol(tillf))
      for(j in 1:nrow(tillf))
      {
        if (!isTRUE(del[j,i]))
          tillf[j,i] <- NA
      }

    tillf <- mice::complete(mice::mice(as.data.frame(tillf),
                                       method = "mean"))
    tillf <- cbind(data[,1], tillf)
    names(tillf) <- names(data)
    return(tillf)
  }
  CooksDistance3M <- function(model)
  {
    cooksd <- stats::cooks.distance(model)
    sample_size <- nrow(data)
#    plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
#    abline(h = (3 * mean(cooksd, na.rm = TRUE)), col="red")  # add cutoff line
 #   text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd > (3 * mean(cooksd, na.rm = TRUE)), names(cooksd),""), col="red")  # add labels
     influential <- as.numeric(names(cooksd)[cooksd > (3 * mean(cooksd, na.rm = TRUE))])
     newData <- data
     newData[influential, ] <- NA
     return(newData)

  }
  RosnerTest <- function()
  {
    newData <- data
    v <-c()

    for (col in newData[-1])
    {
      print(col)

      r <- EnvStats::rosnerTest(col)
      if (r$n.outliers > 0)
      {
        stat <- r$all.stats
        for (index in 1:r$n.outliers)
          v <- c(v, stat[index, "Obs.Num"])
      }
      print(v)
    }
    newData[v,] <- NA
    return(newData)
  }
  data <- data[all.vars(form)]
  ordinaryLM <- lm(form, data)
  dataRosner <- RosnerTest()
  rosner <- lm(form, dataRosner)
  dataCooksLM4N <- CooksDistance4N(ordinaryLM)
  cooksLM4N <- lm(form, dataCooksLM4N)
  dataCooksLM3M <- CooksDistance3M(ordinaryLM)
  cooksLM3M <- lm(form, dataCooksLM3M)
  dataBetas <- DFBetas(ordinaryLM)
  betas <- lm(form,dataBetas )
  dataBetasMiss <- DFBetasMissing(ordinaryLM)
  betasMiss <- lm(form, dataBetasMiss)
  dataChange <- Change()
  changeLM <- lm(form, dataChange)
  dataSD <- SD()
  sdLM <- lm(form, dataSD)
  dataMahalanobis <- Mahalanobis()
  mahalanobisLM <- lm(form, dataMahalanobis)
  dataLeverage <-Leverage(ordinaryLM)
  leverage <- lm(form, dataLeverage)
  robustLM <- robustbase::lmrob(form, data)
  obj <- list(Results = ordinaryLM,
              ordinaryLM = summary(ordinaryLM),
              Rosner = summary(rosner),
              CooksLM4N = summary(cooksLM4N),
              CooksLM3M = summary(cooksLM3M),
              Betas = summary(betas),
              BetasMiss = summary(betasMiss),
              changeLM = summary(changeLM),
              sdLM = summary(sdLM),
              mahalanobisLM = summary(mahalanobisLM),
              Leverage = summary(leverage),
              robustLM = summary(robustLM),
              DataCooksLM4N = dataCooksLM4N,
              DataCooksLM3M = dataCooksLM3M,
              DataBetas = dataBetas,
              DataRosner = dataRosner,
              DataBetasMiss = dataBetasMiss,
              DataChange = dataChange,
              DataSD = dataSD,
              DataMahalanobis = dataMahalanobis,
              DataLeverage = dataLeverage,
              DataCooksLM4N = dataCooksLM4N,
              Data = data)
  class(obj) <- "testOutlier"
  return(obj)

}

#' summaryTestOutlier
#'
#' @param object An object of type testOutliers
#' @param ... Not used now
#' This function prints summary of all lm models with fixed outliers
#' @return nothing
#' @export
summaryTestOutlier <- function(object, ...)
{
  AddN <- function(coeff, N)
  {
     myN <- colnames(coeff)
    coeff <- cbind(coeff, 0)
    coeff[1, ncol(coeff)] <- N
    coeff <- as.data.frame(coeff)
    colnames(coeff) <- c(myN, "N")
    return(coeff)
  }
  print("Ordinary Least Square")
  print(AddN(object$ordinaryLM$coefficients,length(object$ordinaryLM$residuals )))
  print("Ordinary Least Square deleting using Rosner Test")
  print(AddN(object$Rosner$coefficients,length(object$Rosner$residuals )))
  print("Ordinary Least Square deleting using Mahalanobis")
  print(AddN(object$mahalanobisLM$coefficients,length(object$mahalanobisLM$residuals )))
  print("Ordinary Least Square deleting using Leverage - Limit: 2(k+1)/n")
  print(AddN(object$Leverage$coefficients,length(object$Leverage$residuals )))
  print("Ordinary Least Square after removing Cooks distant - Limit: 4 / by N")
  print(AddN(object$CooksLM4N$coefficients,length(object$CooksLM4N$residuals )))
  print("Ordinary Least Square after removing Cooks distant - Limit: 3 times mean of Cooks D")
  print(AddN(object$CooksLM3M$coefficients,length(object$CooksLM3M$residuals )))
  print("Ordinary Least Square after removing any row with DFBETAS - Limit: sqrt(n)")
  print(AddN(object$Betas$coefficients,length(object$Betas$residuals )))
  print("Ordinary Least Square after replacing with mean any cell with DFBETAS - sqrt(n)")
  print(AddN(object$BetasMiss$coefficients,length(object$BetasMiss$residuals )))
  print("Ordinary Least Square after Change/Winsorizing")
  print(AddN(object$changeLM$coefficients,length(object$changeLM$residuals )))
  print("Ordinary Least Square deleteing outside p")
  print(AddN(object$sdLM$coefficients,length(object$sdLM$residuals )))
  print("Robust MM estimation")
  print(AddN(object$robustLM[["coefficients"]], length(object$ordinaryLM$residuals )))
}

#' getDataTestOutlier
#'
#' @param object  An object of type testOutliers
#' @param data A name of the data that should be returned. If data = '?' the function
#' will return all the names in the object that can be used for data (All those
#' names starts with 'DataSomthing')
#'
#' @return a dataframe handled by a specific method
#' @export
getDataTestOutlier <- function(object, data = "?")
{
  if (data == "?")
    print(names(object))
  else
    return(as.data.frame(object[[data]]))
}

#' plotTestOutlier
#' This function plots outlier diagnostics for the model
#' @param x object of type testOutliers
#' @param ... Not used now
#'
#' @return Nothing
#' @export
plotTestOutlier <- function(x, ...)
{
  PlotInfluence <- function(outliers, limit, name)
  {
    plot(outliers, pch="*", cex=2, main="name")  # plot cook's distance
    abline(h = limit, col="red")  # add cutoff line
    text(x=1:length(outliers)+1, y=outliers, labels=ifelse(limit, names(outliers),""), col="red")  # add labels

  }
  print(olsrr::ols_plot_dfbetas(x$Results))
  print(olsrr::ols_plot_cooksd_chart(x$Results))
  print(olsrr::ols_plot_diagnostics(x$Results))
  print(olsrr::ols_test_breusch_pagan(x$Results))
  print(olsrr::ols_test_normality(x$Results))
}

#' Compare datafiles
#'
#' @param object An object of type testOutliers
#' @param what A dataframe to compare
#'
#' @return either names or nothing
#' @export
Compare <- function(object, what = "DataRosner")
{
  UseMethod("Compare", object)
}

#' @export
Compare.testOutlier <- function(object, what="DataRosner")
{
  if (what == "?")
    return(names(object))

  if (what == "All")
  {
  print("Rosner data")
  print(summary(arsenal::comparedf(object$Data, object$DataRosner)))
  print("Mahalanobis")
  print(summary(arsenal::comparedf(object$Data, object$DataMahalanobis)))
  print("Leverage - Limit: 2(k+1)/n")
  print(summary(arsenal::comparedf(object$Data, object$DataLeverage)))
  print("Cooks distant - Limit: 4 / by N")
  print(summary(arsenal::comparedf(object$Data, object$DataCooksLM4N)))
  print("Cooks distant - Limit: 3 times mean of Cooks D")
  print(summary(arsenal::comparedf(object$Data, object$DataCooksLM3M)))
  print("DFBETAS - Limit: sqrt(n)")
  print(summary(arsenal::comparedf(object$Data, object$DataBetas)))
  print("Replacing with mean any cell with DFBETAS - sqrt(n)")
  print(summary(arsenal::comparedf(object$Data, object$DataBetasMiss)))
  print("Change/Winsorizing")
  print(summary(arsenal::comparedf(object$Data, object$DataChange)))
  print("Deleteing outside p")
  print(summary(arsenal::comparedf(object$Data, object$DataSD)))
  }
  else
    print(summary(arsenal::comparedf(object$Data, object[[what]])))


}
