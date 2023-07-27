##
## Prepare data with the Rage lag
##

data <- read.csv("Data/SRdata_5_7_22.csv")

data <- subset(data, !StockKeyLabel %in% c("her.27.1-24a514a", "lin.27.5b"))
data$StockKeyLabel <- as.factor(data$StockKeyLabel)


database  <- list()
res <- list()

for(i in 1:length(levels(data$StockKeyLabel))){
    print(i)
    tmp <- as.data.frame(subset(data, StockKeyLabel == levels(data$StockKeyLabel)[i]))
    ## order by years and then lag
    tmp <- tmp[order(tmp$Year), ]
    rlag <- unique(tmp$RecruitmentAge)
    if(rlag == 0){
        tmp <- data.frame( StockKeyLabel = tmp$StockKeyLabel, AssessmentYear = tmp$AssessmentYear,
                          Year = tmp$Year, SSB = tmp$SSB, Recruitment = tmp$Recruitment, F = tmp$F)
    }else{
        tmp <- data.frame(StockKeyLabel = tmp$StockKeyLabel, AssessmentYear = tmp$AssessmentYear,
                          Year = tmp$Year, SSB = tmp$SSB, Recruitment = c(tmp$Recruitment[-(1:rlag)],
                                                                          rep(NA, rlag)), F = tmp$F)
    }
    ##
    #tmp <- subset(tmp, !is.na(tmp$Recruitment))
    if(nrow(tmp) <= 20){
        NA   
    }else{
    ##
    database <- rbind(database, tmp)
    }
    }

write.csv(database, "Data/SR_CSEdata.csv")
