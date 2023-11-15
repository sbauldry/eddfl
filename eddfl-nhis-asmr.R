### Purpose: Prepare NHIS-LMF for estimates of age-specific mortality rates
### Author:  S Bauldry
### Date:    September 29, 2022

### Loading packages
setwd("~/desktop")
library(ipumsr)
library(tidyverse)
library(lubridate)
library(clock)

### Read data extract
ddi <- read_ipums_ddi("nhis_00003.xml")
d1 <- read_ipums_micro(ddi)

### Prepare variables for analysis
d1 <- d1 %>%
  rename_with(tolower) %>%
  filter(year >= 2000, age >= 65, mortelig == 1) %>%
  mutate(
    # prepare covariates
    fem = ifelse(sex == 2, 1, 0), 
    edu = case_when(
      educ <= 116 ~ 1,  #less than high school
      educ <= 303 ~ 2,  # high school, GED, some college, vocational training, or AA degree
      educ <= 503 ~ 3), # BA or higher
    
    # set top-coded age at 85 to missing
    age = ifelse(age < 85, age, NA),
    
    # set missing birth month to mid-year
    bmo = ifelse(birthmo < 97, birthmo, 7),
    
    # set missing values for birth year
    byr = ifelse(birthyr < 9997, birthyr, NA)
 )

### Set bottom-coded birth year to missing
for(i in 1:16) {
  bry = 1914 + i
  sry = 1998 + i
  d1$byr[d1$year == sry & d1$byr == bry] = NA
}

### Calculate dates and age on Jan 1
d1 <- d1 %>%
  mutate(
    
    # interview date
    intyr = ifelse(intervwyr == 9998, year, intervwyr),
    doi = date_build(intyr, intervwmo, 15),
    
    # birth date
    byr = ifelse(is.na(byr), intyr - age, byr),
    dob = date_build(byr, bmo, 15),
    
    # date of death or censoring
    dth = ifelse(mortstat == 1, 1, 0),
    dod = case_when(
      dth == 0 ~ as_date( date_build(2015, 12, 31) ),
      mortdodq == 1 ~ as_date( date_build(mortdody, 2, 15) ),
      mortdodq == 2 ~ as_date( date_build(mortdody, 5, 15) ),
      mortdodq == 3 ~ as_date( date_build(mortdody, 8, 15) ),
      mortdodq == 4 ~ as_date( date_build(mortdody, 11, 15) )),
    
    # fix 121 cases with interview date after death date
    dod = case_when(
      doi <= dod ~ dod,
      doi >  dod ~ doi + 1),
    
    # age on Jan 1
    jage = as.numeric( round( (date_build(intyr, 1, 1) - dob)/365.25 ) )
 )

### Select analysis sample and variables before expanding
# baseline sample size
d1 <- d1 %>% filter(jage >= 65)
length(d1$nhispid)

# drop weights of 0
d1 <- d1 %>% filter(mortwt != 0)
length(d1$nhispid)

# drop missing date of death
d1 <- d1 %>% drop_na(dod)
length(d1$nhispid)

# drop missing demographic variables
d1 <- d1 %>% drop_na(c(fem, edu))
length(d1$nhispid)

# select analysis variables
d1 <- d1 %>% select(nhispid, intyr, fem, edu, jage, doi, dod, dth, mortwt, 
                    strata, psu)


### Expand to person-year format
d2 <- d1 %>%
  mutate(ny = ifelse(intyr < year(dod), year(dod) + 1 - intyr, 1)) %>%
  uncount(ny)
length(d2$nhispid)


### Prepare person-year data for survival models
d3 <- d2 %>%
  group_by(nhispid) %>%
  mutate(
    
    # calculate current year
    cyear = row_number() + intyr - 1,
    
    # calculate age in current year
    cage = jage + row_number() - 1,
    cage = ifelse(cage > 85, 85, cage),
    
    # create indicator for died in current year
     died = ifelse(dth == 1 & cyear == max(cyear), 1, 0),
    
    # create exposure variable in current year
    cyrb = date_build(cyear, 1, 1),
    cyre = date_build(cyear, 12, 31),
    cerd = ifelse(cyear == year(dod) & died == 1, (dod - cyrb + 1)/365.25, 1),
    cerd = ifelse(cyear == intyr & died == 0, (cyre - doi + 1)/365.25, cerd)
  ) %>%
  select(nhispid, fem, edu, cage, cerd, died, mortwt, strata, psu)

### Save data for bootstrapping
write_csv(d3, "eddfl-asmr-data.csv")
