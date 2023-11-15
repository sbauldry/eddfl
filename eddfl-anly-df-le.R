### Purpose: Conduct analysis of education and 2FLEs
### Author:  S Bauldry
### Date:    Nov 7, 2023

### Load packages and set working directory
setwd("~/desktop")
library(tidyverse)
library(ggpubr)
library(survey)
library(weights)


### Read prepared age-specific dual-function and mortality data
### Combine strata for 9 cases
hrs  <- read_csv("eddfl-asdf-data.csv", col_types = list(agei = "f", edu = "f")) |>
  mutate(strt = ifelse(raestrat > 60, 60, raestrat))
nhis <- read_csv("eddfl-asmr-data.csv") 


### Proportion of observations based on proxy reports in HRS
prop.table( table(hrs$pxy) )


### Descriptive statistics (Table 1)
# HRS
length(hrs$hhidpn)
n_distinct(hrs$hhidpn)

wpct(hrs$fem, weight = hrs$wtcrnh)
wpct(hrs$edu, weight = hrs$wtcrnh)
wtd.mean(hrs$age, weight = hrs$wtcrnh)
sqrt(wtd.var(hrs$age, weight = hrs$wtcrnh))

# NHIS
nhis_pl <- nhis %>%
  group_by(nhispid) %>%
  filter(row_number() > (n() - 1))

length(nhis_pl$nhispid)

wpct(nhis_pl$fem, weight = nhis_pl$mortwt)
wpct(nhis_pl$edu, weight = nhis_pl$mortwt)
wtd.mean(nhis_pl$cage, weight = nhis_pl$mortwt)
sqrt(wtd.var(nhis_pl$cage, weight = nhis_pl$mortwt))



### Plot age-group specific dual-function rates by gender and education (Figure 1)
dfwp <- function(d1) {
  des <- svydesign(id = ~raehsamp, strata = ~strt, weights = ~wtcrnh, nest = T, data = d1)
  df <- data.frame( expand.grid( agei = seq(65, 85, 5), edu = c(1, 2, 3) ) )
  df$agei <- factor(df$agei)
  df$edu  <- factor(df$edu) 
  m1  <- svyglm(df1 ~ agei + edu + agei*edu, data = d1, family = quasibinomial, design = des)
  pdf1 <- predict(m1, newdata = df, type = "response")
  pdf <- cbind(df, pdf1)
}

# dual-function rates for women
dfp_fem <- hrs %>%
  filter(fem == 1) %>%
  dfwp() %>%
  rename(fdfp = response, fdfse = SE) %>%
  mutate(pdf = fdfp,
         plb = (fdfp - 1.96*fdfse),
         pub = (fdfp + 1.96*fdfse),
         agei = as.numeric(agei))

# dual-function rates for men
dfp_mal <- hrs %>%
  filter(fem == 0) %>%
  dfwp() %>%
  rename(mdfp = response, mdfse = SE) %>%
  mutate(pdf = mdfp,
         plb = (mdfp - 1.96*mdfse),
         pub = (mdfp + 1.96*mdfse),
         agei = as.numeric(agei))

# generate figure
figa <- function(dfn, tit, df) {
  fig <- ggplot(data = dfn, mapping = aes(x = agei, y = 100*pdf, color = edu)) +
    geom_line() +
    geom_point(size = 2) +
    geom_ribbon(mapping = aes(ymin = 100*plb, ymax = 100*pub), alpha = 0.2, linetype = 0) +
    scale_color_manual(values = c("1" = "blue", "2" = "dark green", "3" = "dark red")) +
    scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 100), name = "percent dual functional") +
    scale_x_continuous(breaks = 1:5, name = "age interval", labels = c("65-69", "70-74", "75-79", "80-84", "85+")) +
    labs(title = tit) +
    guides(fill = "none", color = "none") +
    theme_light(base_size = 15) 
  return(fig)
}

fig1a <- figa(dfp_fem, "Panel A: Women")
fig1b <- figa(dfp_mal, "Panel B: Men")
fig1 <- ggarrange(fig1a, fig1b) 
ggsave("eddfl-fig1.jpg", plot = fig1, width = 10, height = 7, units = "in")



### Differences across education by gender and age-group (Table 2)
dfp_fem_dif <- dfp_fem |>
  select(agei, edu, fdfp, fdfse) |>
  pivot_wider(names_from = edu, values_from = c(fdfp, fdfse)) |>
  mutate(d1 = fdfp_2 - fdfp_1,
         d1se = sqrt( fdfse_2^2 + fdfse_1^2 ),
         d1pv = 2*(1 - pnorm( abs(d1/d1se) ) ),
         d2 = fdfp_3 - fdfp_1,
         d2se = sqrt( fdfse_3^2 + fdfse_1^2 ),
         d2pv = 2*(1 - pnorm( abs(d2/d2se) ) )) |>
  select(d1, d1pv, d2, d2pv)
dfp_fem_dif

dfp_mal_dif <- dfp_mal |>
  select(agei, edu, mdfp, mdfse) |>
  pivot_wider(names_from = edu, values_from = c(mdfp, mdfse)) |>
  mutate(d1 = mdfp_2 - mdfp_1,
         d1se = sqrt( mdfse_2^2 + mdfse_1^2 ),
         d1pv = 2*(1 - pnorm( abs(d1/d1se) ) ),
         d2 = mdfp_3 - mdfp_1,
         d2se = sqrt( mdfse_3^2 + mdfse_1^2 ),
         d2pv = 2*(1 - pnorm( abs(d2/d2se) ) )) |>
  select(d1, d1pv, d2, d2pv)
dfp_mal_dif




### Age-65 life expectancy and dual-function life expectancy (Table 3)
# function to read bootstrap data and extract results
bsres <- function(fn) {
  bsd <- read_csv(fn)
  est <- c( mean(bsd$e65) + 65, sd(bsd$e65), 
            mean(bsd$dfe65) + 65, sd(bsd$dfe65), 
            mean(bsd$dfe65)/mean(bsd$e65) )
  return(est)
}

le_lh_fem <- bsres("eddfl-age65-le-eg-11.csv")
le_hs_fem <- bsres("eddfl-age65-le-eg-21.csv")
le_ba_fem <- bsres("eddfl-age65-le-eg-31.csv")

le_lh_mal <- bsres("eddfl-age65-le-eg-10.csv")
le_hs_mal <- bsres("eddfl-age65-le-eg-20.csv")
le_ba_mal <- bsres("eddfl-age65-le-eg-30.csv")

le_fem <- rbind( c(1, le_lh_fem), c(2, le_hs_fem), c(3, le_ba_fem) )
le_mal <- rbind( c(1, le_lh_mal), c(2, le_hs_mal), c(3, le_ba_mal) )
colnames(le_fem) <- c("edu", "le", "lese", "dfle", "dflese", "pdf")
colnames(le_mal) <- c("edu", "le", "lese", "dfle", "dflese", "pdf")
le_fem
le_mal

# calculate differences
le_fem_dif <- as_tibble(le_fem) |>
  select(-pdf) |>
  pivot_wider(names_from = edu, values_from = c(le, lese, dfle, dflese)) |>
  mutate(d1 = le_2 - le_1,
         d1se = sqrt( lese_2^2 + lese_1^2 ),
         d1pv = 2*(1 - pnorm( abs(d1/d1se) ) ),
         d2 = le_3 - le_1,
         d2se = sqrt( lese_3^2 + lese_1^2 ),
         d2pv = 2*(1 - pnorm( abs(d2/d2se) ) ),
         d3 = dfle_2 - dfle_1,
         d3se = sqrt( dflese_2^2 + dflese_1^2 ),
         d3pv = 2*(1 - pnorm( abs(d3/d3se) ) ),
         d4 = dfle_3 - dfle_1,
         d4se = sqrt( dflese_3^2 + dflese_1^2 ),
         d4pv = 2*(1 - pnorm( abs(d4/d4se) ) )) |>
  select(d1, d1pv, d2, d2pv, d3, d3pv, d4, d4pv)
le_fem_dif

le_mal_dif <- as_tibble(le_mal) |>
  select(-pdf) |>
  pivot_wider(names_from = edu, values_from = c(le, lese, dfle, dflese)) |>
  mutate(d1 = le_2 - le_1,
         d1se = sqrt( lese_2^2 + lese_1^2 ),
         d1pv = 2*(1 - pnorm( abs(d1/d1se) ) ),
         d2 = le_3 - le_1,
         d2se = sqrt( lese_3^2 + lese_1^2 ),
         d2pv = 2*(1 - pnorm( abs(d2/d2se) ) ),
         d3 = dfle_2 - dfle_1,
         d3se = sqrt( dflese_2^2 + dflese_1^2 ),
         d3pv = 2*(1 - pnorm( abs(d3/d3se) ) ),
         d4 = dfle_3 - dfle_1,
         d4se = sqrt( dflese_3^2 + dflese_1^2 ),
         d4pv = 2*(1 - pnorm( abs(d4/d4se) ) )) |>
  select(d1, d1pv, d2, d2pv, d3, d3pv, d4, d4pv)
le_mal_dif

  


