### Purpose: Prepare bootstrap age-50 life expectancy estimates
### Author:  S Bauldry
### Date:    Nov 7, 2023

### load packages and set working directory
setwd("~/desktop")
library(tidyverse)
library(survival)
library(survey)


### load prepared HRS and NHIS data
nhis <- read_csv("eddfl-asmr-data.csv") 
hrs  <- read_csv("eddfl-asdf-data.csv") 


### Function for life table calculations
lifetable <- function(x) {
  
  # calculate qx and lx
  for(i in 1:nrow(x)) {
    if(i < nrow(x)) {
      x[i,"qx"] = (x[i,"nx"]*x[i,"mx"])/(1 + (x[i,"nx"] - x[i,"ax"])*x[i,"mx"])
    }
    if(i == nrow(x)) {
      x[i,"qx"] = 1
    }
    
    if(i == 1) {
      x[i,"lx"] = 100000
    }
    if(i > 1) {
      x[i,"lx"] = x[i-1,"lx"]*(1 - x[i-1,"qx"])
    }
  }
  
  # calculate dx, Lx, and dLx
  for(i in 1:nrow(x)) {
    if(i < nrow(x)) {
      x[i,"dx"] = x[i,"lx"] - x[i+1,"lx"]
      x[i,"Lx"] = x[i,"nx"]*x[i+1,"lx"] + x[i,"dx"]*x[i,"ax"]
    }
    if(i == nrow(x)) {
      x[i,"dx"] = x[i,"lx"]
      x[i,"Lx"] = x[i,"lx"]/x[i,"mx"]
    }
    x[i,"dLx"] = x[i,"Lx"]*(x[i,"dfp"])
  }
  
  # calculate Tx and dTx
  for(i in nrow(x):1) {
    if(i < nrow(x)) {
      x[i,"Tx"]  = x[i+1,"Tx"] + x[i,"Lx"]
      x[i,"dTx"] = x[i+1,"dTx"] + x[i,"dLx"]
    } 
    if(i == nrow(x)) {
      x[i,"Tx"]  = x[i,"Lx"]
      x[i,"dTx"] = x[i,"dLx"]
    }
  }
  
  # calculate ex and dex
  for(i in 1:nrow(x)) {
    x[i,"ex"]  = x[i,"Tx"]/x[i,"lx"]
    x[i,"dex"] = x[i,"dTx"]/x[i,"lx"]
  }
  
  # return age-50 ex and dex
  e <- x[1, c("agei", "ex", "dex")]
  return(e)
}


### Function to calculate age-50 life expectancy for bootstrapped sample
le_est <- function(d1, d2) {
  
  # obtain estimates of age-specific mortality rates
  nuid <- unique(d1$nhispid)
  nsam <- tibble( sample(nuid, size = length(nuid), replace = T) )
  colnames(nsam) <- "nhispid"
  nbs <- nsam %>% left_join(d1, by = "nhispid", relationship = "many-to-many")
  nhis_des <- svydesign(id = ~1, weights = ~mortwt, data = nbs)
  nm  <- svysurvreg(Surv(cerd, died) ~ cage, dist = "exponential", data = nbs, design = nhis_des)
  aslp <- predict(nm, data.frame(cage = 65:94), type = "lp")
  ashr <- exp(-aslp)
  ai <- c( sort( rep(1:4, 5) ), rep(5, 10) )
  mx <- tapply(ashr, ai, mean)
  
  # obtain estimates of age-specific dual function rates
  huid <- unique(d2$hhidpn)
  hsam <- tibble( sample(huid, size = length(huid), replace = T) )
  colnames(hsam) <- "hhidpn"
  hbs <- hsam %>% left_join(d2, by = "hhidpn", relationship = "many-to-many")
  hrs_des <- svydesign(id = ~1, weights = ~wtcrnh, data = hbs)
  m1 <- svyglm(df1 ~ as.factor(agei), data = hbs, family = quasibinomial, design = hrs_des)
  aged <- data.frame( agei = seq(65, 85, 5) )
  asdf <- cbind( aged, predict(m1, newdata = aged, type = "response"))

  # calculate age-50 total and dual-function life expectancy
  agei <- seq(65, 85, 5)
  nx   <- c( rep(5, 4), 10 )
  ax   <- c( rep(2.5, 4), 5 )
  dfp <- asdf$response  
  lt  <- tibble(agei, nx, ax, mx, dfp)
  e1  <- lifetable(lt)
  
  # combine estimates
  ec <- c( e1[[1,1]], e1[[1,2]], e1[[1,3]] )
    
  # return estimates
  return(ec)
}



### Set seed for reproducing bootstrap results and number of bootstrap samples
set.seed(585399)
nb <- 500

### Bootstrap estimates by education and gender
for(s in 1:3) {
  for(g in 0:1) {
    nhis_sub <- nhis %>% filter(edu == s & fem == g)
    hrs_sub  <- hrs  %>% filter(edu == s & fem == g)
    est <- matrix(NA, nb, 3)
    for(i in 1:nb) {
      print( c(s, g, i) )
      be <- le_est(nhis_sub, hrs_sub)
      est[i,] <- c(i, be[2], be[3])
    }
    colnames(est) <- c("bsam", "e65", "dfe65")
    est <- data.frame(est)
    fn  <- paste("eddfl-age65-le-eg-", paste0(s,g) , ".csv", sep = "") 
    write_csv(est, fn)
  }
}
