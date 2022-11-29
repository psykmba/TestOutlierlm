library(devtools)
use_r("TestOutliers")


library(readr)
devtools::install_github("psykmba/DoPsychometric", force = T)
library(DoPsychometric)


KallesTest200625Scales <- read_csv("~/Dropbox/Forskning/SPSS filer/Personlighetstest/Kalle/KallesTest200625Scales.txt")
PWBScales <- read_csv("~/Dropbox/Forskning/SPSS filer/Personlighetstest/Kalle/PWBScales.txt")
View(PWBScales)
testData <- dplyr::inner_join(KallesTest200625Scales[c(1,5,6, 16:37)], PWBScales[c(1,16:21)], by = c("ID"))
View(testData)
names(testData)
pObject <- DoPsychometric::GetPsychometric(testData, scaleNames = c("E", "A", "C", "N", "O"),
                                           itemList = list(c(4, 10, 22, 12),
                                                           c(8, 20, 21, 24),
                                                           c(6, 13, 14, 25),
                                                           c(5, 9, 15, 16),
                                                           c(7, 11, 17, 23)),
                                           responseScale = list(c(0,64)),
                                           reverse = F,
                                           itemLength = 1)
testData <- DoPsychometric::getData(pObject)

usethis::use_data(testData)

res <- testOutliers(SelfAcceptance ~ E + C + A + N + O, testData)
summary.testOutlier(res)
testData[1, "PurposeInLife"]= 80
getData.testOutlier(res, "DataBetas")
res2 <- testOutliers(SelfAcceptance ~ E + C + A + N + O, getData.testOutlier(res, "DataBetas"), limit = .01)
outliers::grubbs.test(testData$E)
