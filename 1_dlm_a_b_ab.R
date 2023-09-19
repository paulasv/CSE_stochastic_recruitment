###############
## Paula
## univariate multistock, dlm, cpp
###############


library(dlm)
library(mgcv)
library(ggplot2); theme_set(theme_classic())
library(reshape)
library(ggpubr)
library(ggrepel)

library(dplyr)
library(forcats)

## data
dat <- read.csv("Data/SR_CSEdata_c.csv")
names(dat)[names(dat) == "StockKeyLabel"] <- "ID"
names(dat)[names(dat) == "Recruitment"] <- "Recruits"

## loop
ll_res <- NULL
fit_res <- NULL

stocks <- unique(dat$ID)

for(i in 1:length(stocks)){
        stock_id <- stocks[i]
        print(stocks[i])
        sub_dat <- subset(dat, ID == stock_id)
        sub_dat <- sub_dat[!is.na(sub_dat$Recruits),]
        ## make sure data ordered correctly
        sub_dat <- sub_dat[order(sub_dat$Year),]
        n <- nrow(sub_dat)
        sub_dat$y <- with(sub_dat, log(Recruits/SSB))
        years <- sub_dat$Year
        ##-----
        ## DLM 
        ##-----
        ## time-varying a
        build_tv <- function(u) {
            mod <- dlmModReg(sub_dat$SSB, dV = exp(u[1]))
            diag(mod$W) <- c(exp(u[2]), 0)
            return(mod)
        }
        fit_tv <- dlmMLE(sub_dat$y, parm = rep(log(0.2), 2),
                         build_tv, method = "L-BFGS-B")
        print("time varying a")
        print(fit_tv$convergence)
        #fit_tv <- dlmMLE(sub_dat$y, parm = rep(log(0.2), 2),
                         #build_tv)
        ## time-invariant
        build_ti <- function(u) {
            mod <- dlmModReg(sub_dat$SSB, dV = exp(u[1]))
            diag(mod$W) <- c(0, 0)
            return(mod)
        }
        fit_ti <- dlmMLE(sub_dat$y, parm = log(0.2),
                         build_ti, method = "L-BFGS-B")
        print("time invariant")
        print(fit_ti$convergence)
        #fit_ti <- dlmMLE(sub_dat$y, parm = log(0.2),
                                        #build_ti)
        dlmti <- dlmSmooth(sub_dat$y, build_ti(fit_ti$par))$s[1,1]
        dlmtib <- dlmSmooth(sub_dat$y, build_ti(fit_ti$par))$s[1,2]
        ## time-varying b
        build_tvb <- function(u) {
            mod <- dlmModReg(sub_dat$SSB, dV = exp(u[1]))
            diag(mod$W) <- c(0, exp(u[2]))
            return(mod)
        }
        
        fit_tv.b <- dlmMLE(sub_dat$y, parm = rep(log(0.2), 2),
                           build_tvb, method = "L-BFGS-B")
        print("time varying b")
        print(fit_tv.b$convergence)
        ## time covaring a b
        mod <- dlmModReg(sub_dat$SSB)
        build_tvab <- function(theta){
            mod1 <- mod
            L <- matrix(c( exp(theta[1]),  exp(theta[2]), 0,  exp(theta[3])), nrow = 2)
            mod1$W <-L %*% t(L)
            mod1$V[1, 1] <- exp(theta[4])
            return(mod1)
        }

        fit_tv.ab <- dlmMLE(sub_dat$y, parm = rep(log(0.2), 4),
                            build_tvab, method = "L-BFGS-B")
        print("time varying ab")
        print(fit_tv.ab$convergence)

        ##
        if((fit_ti$convergence == 0) & (fit_tv$convergence == 0)){
            llti <- fit_ti$value # dlm calculates the negative loglikelihood
            lltv <- fit_tv$value
            lltv.b <- fit_tv.b$value
            lltv.ab <- fit_tv.ab$value
            
            D <- 2*(llti - lltv)
            pvalue <- 1-pchisq(D, 1)

            D.b <- 2*(llti - lltv.b)
            pvalue.b <- 1-pchisq(D.b, 1)

            D.ab <- 2*(llti - lltv.ab)
            pvalue.ab <- 1-pchisq(D.ab, 1)
           
            AICti <- 2 * llti + 2 * 3 #*(n/(n-3-1))
            AICtv <- 2 * lltv + 2 * 4 #*(n/(n-4-1))
            AICtv.b <- 2 * lltv.b + 2 * 4 #*(n/(n-4-1))
            AICtv.ab <- 2 * lltv.ab + 2 * 6 #*(n/(n-6-1))
            ##
            dlmres <- data.frame(stock = stock_id,
                                 method = "dlm",
                                 llti = llti,
                                 lltv = lltv,
                                 lltv.b = lltv.b,
                                 lltv.ab = lltv.ab,
                                 D = D,
                                 df = 1,
                                 pvalue = pvalue,
                                 pvalue.b = pvalue.b,
                                 pvalue.ab = pvalue.ab,
                                 AICti = AICti,
                                 AICtv = AICtv,
                                 AICtv.b = AICtv.b,
                                 AICtv.ab = AICtv.ab
                                 )
            ##
            smooth <- dlmSmooth(sub_dat$y, build_tv(fit_tv$par))
            se_dlm <- unlist(lapply(dlmSvd2var(smooth$U.S, smooth$D.S), function(x){sqrt(x[1,1])}))
            adlm <- data.frame(ahat = dropFirst(smooth$s[,1]),
                               ase = dropFirst(se_dlm))
            se_dlmb <- unlist(lapply(dlmSvd2var(smooth$U.S, smooth$D.S), function(x){sqrt(x[2,2])}))
            bdlm <- data.frame(bhat = dropFirst(smooth$s[,2]),
                               bse = dropFirst(se_dlmb))
            ##
            smooth.b <- dlmSmooth(sub_dat$y, build_tvb(fit_tv.b$par))
            se_dlm.b <- unlist(lapply(dlmSvd2var(smooth.b$U.S, smooth.b$D.S), function(x){sqrt(x[1,1])}))
            adlm.b <- data.frame(ahat = dropFirst(smooth.b$s[,1]),
                                 ase = dropFirst(se_dlm.b))
            se_dlmb.b <- unlist(lapply(dlmSvd2var(smooth.b$U.S, smooth.b$D.S), function(x){sqrt(x[2,2])}))
            bdlm.b <- data.frame(bhat = dropFirst(smooth.b$s[,2]),
                               bse = dropFirst(se_dlmb.b))
             ##
            smooth.ab <- dlmSmooth(sub_dat$y, build_tvab(fit_tv.ab$par))
            se_dlm.ab <- unlist(lapply(dlmSvd2var(smooth.ab$U.S, smooth.ab$D.S), function(x){sqrt(x[1,1])}))
            adlm.ab <- data.frame(ahat = dropFirst(smooth.ab$s[,1]),
                                  ase = dropFirst(se_dlm.ab))
            se_dlmb.ab <- unlist(lapply(dlmSvd2var(smooth.ab$U.S, smooth.ab$D.S), function(x){sqrt(x[2,2])}))
            bdlm.ab <- data.frame(bhat = dropFirst(smooth.ab$s[,2]),
                               bse = dropFirst(se_dlmb.ab))
            ##
            names(adlm) <- c("ahat", "ase")
        }else{
            dlmres <- data.frame(stock = stock_id,
                                 method = "dlm",
                                 llti = NA,
                                 lltv = NA,
                                 lltv.b = NA,
                                 lltv.ab = NA,
                                 D = NA,
                                 df = NA,
                                 pvalue = NA,
                                 pvalue.b = NA,
                                 pvalue.ab = NA,
                                 AICti = NA,
                                 AICtv = NA,
                                 AICtv.b = NA,
                                 AICtv.ab = NA)
           adlm <- data.frame(ahat = rep(NA, n), ase = rep(NA, n))
        }
        ##
        ll_res <- rbind(ll_res,  dlmres)
        
        res <- cbind(sub_dat, dlmti, dlmtib, adlm,bdlm, "bt" = adlm.b, "bt" = bdlm.b, "abt" = adlm.ab, "abt" = bdlm.ab)
        fit_res <- rbind(fit_res, res)
    }


# time-invariant
idx <- with(ll_res, which(AICti<AICtv & AICti<AICtv.b & AICti< AICtv.ab))
ti_res <- ll_res[idx,]
# time-varing a
idx <- with(ll_res, which(AICtv<AICti & AICtv<AICtv.b & AICtv< AICtv.ab))
tv_res <- ll_res[idx,]
# time-varing b
idx <- with(ll_res, which(AICtv.b<AICti & AICtv.b<AICtv & AICtv.b< AICtv.ab))
tv_res.b <- ll_res[idx,]
# time-varing ab
idx <- with(ll_res, which(AICtv.ab<AICti & AICtv.ab<AICtv & AICtv.ab< AICtv.b))
tv_res.ab <- ll_res[idx,]

##-----------------------------
## FIGURES
##-----------------------------

data_tva <- subset(fit_res, ID %in% tv_res$stock)

data_ti <- cbind(aggregate(data = data_tva, ahat~ID, FUN = mean), dlmti = aggregate(data = data_tva, dlmti~ID, FUN = mean)[,2])

## Figure 1
data_sub <- subset(data_tva, Year>1955)
stocks <- read.csv("Data/stock_info.csv")
stocks$ID <- stocks$stocks
data.cat <- data_sub %>%
    left_join(stocks)

data.cat <- data.cat%>%
    mutate(Region = fct_relevel(Region, c("Northern CSE", "Central CSE", "Southern CSE", "NEA wide area")))%>%
    mutate(ID = fct_reorder(ID,  as.numeric(Region)))

data.cat <- data.cat%>%
    mutate(ID = fct_relevel(ID, c("aru.27.5b6a", "had.27.46a20", "whg.27.6a", "her.27.6a7bc", "had.27.6b", "cod.27.7a", "whg.27.7a", "sol.27.7a", "her.27.nirs", "cod.27.7e-k", "whg.27.7b-ce-k", "sol.27.7fg", "meg.27.7b-k8abd", "mon.27.78abd", "reg.27.561214", "whb.27.1-91214", "dgs.27.nea", "mac.27.nea")))

jpeg("Plots/Figure1.jpg", width = 4200, height = 4900, res = 420)
ggplot(data.cat, aes(Year, ahat))+
    geom_line(colour = "steelblue3")+
    geom_ribbon(aes(ymin = ahat - ase, ymax = ahat + ase), fill = "steelblue1", alpha = .2)+
    #geom_hline(yintercept = 0, colour = "gray", linetype = 3)+
    geom_line(aes(Year, dlmti), colour = "chartreuse2", linetype = 2)+
    scale_x_continuous(breaks=seq(1960,2020,10))+
    facet_wrap(Region~ID, ncol = 3,  scales = "free_y")+
    ylab(expression(a[t]))+
    labs(colour = "")
dev.off()

## Figure 2
ts <- ggplot(data_ti, aes(ahat, dlmti))+
    geom_point()+
    geom_abline(slope = 1, linetype = 3)+
    geom_text_repel(aes(label = ID))+
    ylim(-0.01,NA)+xlim(-0.01,NA)+
    xlab(expression(mean(a[t])))+ ylab("maximum-productivity parameter (a)")+ ggtitle("(A)")+
    theme(axis.title = element_text(size=17), axis.text = element_text(size = 12))

data_tva_2018 <- subset(data_tva, Year>2016)
data_ti <- cbind(aggregate(data = data_tva_2018, ahat~ID, FUN = mean), dlmti = aggregate(data = data_tva_2018, dlmti~ID, FUN = mean)[,2])
ts_16 <- ggplot(data_ti, aes(ahat, dlmti))+
    geom_point()+
    geom_abline(slope = 1, linetype = 3)+
    #geom_text(aes(label = ID), vjust = -.7, size = 3.2)+
    geom_text_repel(aes(label = ID))+
    ylim(-0.01,NA)+xlim(-0.3,NA)+
    xlab(expression(mean(a[t])))+ ylab("maximum-productivity parameter (a)")+ ggtitle("(B)")+
    theme(axis.title = element_text(size=17), axis.text = element_text(size = 12))


jpeg("Plots/Figure2.jpg", width = 4000, height = 1800, res = 300)
ggarrange(ts,ts_16)
dev.off()



## Figure 3
data_tvb <- subset(fit_res, ID %in% tv_res.b$stock)

jpeg("Plots/Figure3.jpg", width = 2900, height = 1900, res = 420)
ggplot(data_tvb, aes(Year, bt.bhat))+
    geom_line(colour = "tomato")+
    geom_ribbon(aes(ymin = bt.bhat - bt.bse, ymax = bt.bhat + bt.bse), fill = "tomato", alpha = .2)+
    geom_line(aes(Year, dlmtib), colour = "orange", linetype = 2)+
    scale_x_continuous(breaks=seq(1960,2020,10))+
    facet_wrap(vars(ID), ncol = 2, scales = "free_y")+
    ylab(expression(b[t]))
dev.off()


##############
## Suplementary plots
###########

data <- subset(fit_res, Year>1955)


jpeg("Plots/plot_Sup1.jpg", width = 5200, height = 5600, res = 420)
ggplot(data, aes(Year, ahat))+
    geom_line(colour = "steelblue3")+
    geom_ribbon(aes(ymin = ahat - ase, ymax = ahat + ase), fill = "steelblue1", alpha = .2)+
    geom_line(aes(Year, dlmti), colour = "black", linetype = 3)+
    geom_line(aes(Year, bt.ahat), linetype = 5, colour = "tomato3")+
    geom_ribbon(aes(ymin = bt.ahat - bt.ase, ymax = bt.ahat + bt.ase), fill = "tomato1", alpha = .2)+
    geom_line(aes(Year, abt.ahat), linetype = 5, colour = "darkgoldenrod1")+
    geom_ribbon(aes(ymin = abt.ahat - abt.ase, ymax = abt.ahat + abt.ase), fill = "darkgoldenrod1", alpha = .2)+
    scale_x_continuous(breaks=seq(1960,2020,10))+
    facet_wrap(vars(ID), ncol = 4, scales = "free_y")+
    ylab(expression(a[t]))
dev.off()


jpeg("Plots/plot_Sup2.jpg", width = 5200, height = 5600, res = 420)
ggplot(data, aes(Year, bhat))+
    geom_line(colour = "steelblue3")+
    geom_ribbon(aes(ymin = bhat - bse, ymax = bhat + bse), fill = "steelblue1", alpha = .2)+
    geom_line(aes(Year, dlmtib), colour = "black", linetype = 3)+
    geom_line(aes(Year, bt.bhat), linetype = 5, colour = "tomato3")+
    geom_ribbon(aes(ymin = bt.bhat - bt.bse, ymax = bt.bhat + bt.bse), fill = "tomato1", alpha = .2)+
    geom_line(aes(Year, abt.bhat), linetype = 5, colour = "darkgoldenrod1")+
    geom_ribbon(aes(ymin = abt.bhat - abt.bse, ymax = abt.bhat + abt.bse), fill = "darkgoldenrod1", alpha = .2)+
    scale_x_continuous(breaks=seq(1960,2020,10))+
    facet_wrap(vars(ID), ncol = 4, scales = "free_y")+
    ylab(expression(b[t]))
dev.off()
