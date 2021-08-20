########################################
##  TITLE OF SCRIPT GOES HERE  ##
########################################
# 
# AUTHOR: SVD
# 
# EMAIL: svd14@case.edu
# 
# DATE: 2021-08-20
# 
# SCRIPT PURPOSE: Random effects metanalysis using logHR
#   
# SCRIPT OTHER DETAILS: random effects meta-analysis using logHR 
# to pool ES at 12, 24, 36 and 48 months.


# SET WORKING DIRECTORY ---

# IMPORT PACKAGE LIBRARIES NEEDED ---

cat("IMPORTING THESE PACKAGES... \n\n", sep = "")
packages <- c("tidyverse","survival","metafor","ckbplotr",
              "readxl")
n_packages <- length(packages)

# install missing packages 

new.packages <- packages[!(packages %in% installed.packages())]

if(length(new.packages)){
  install.packages(new.packages)
}

# load all libraries needed

for(n in 1:n_packages){
  cat("Loading Library #", n, "of", n_packages, "...currently loading: ", packages[n], "\n", sep = "")
  lib.load <- paste("library(\"", packages[n], "\")", sep = "")
  eval(parse(text = lib_load))
}

# SETTING OPTIONS ---

cat("SETTING OPIONS... \n\n", sep = "")
options(scipen = 999)
options(encoding = "UTF-8")

#---------------------------------------------------
# START CODING HERE
#---------------------------------------------------


# get the dataset for analysis - 


# getting for 12 months - 

df12 = readxl::read_excel('F:\\GLP1_agonists\\analysis\\cox_models\\data.xlsx',
                        sheet = "12months")

glimpse(df12)

res12 = metafor::rma.uni(yi = logHR, sei = se, data = df12,
                         method = 'DL')

summary(res12)

plot(res12, exp = T)

# getting for 24 months - 

df24 = readxl::read_excel('F:\\GLP1_agonists\\analysis\\cox_models\\data.xlsx',
                          sheet = "24months")

glimpse(df24)

res24 = metafor::rma.uni(yi = logHR, sei = se, data = df24,
                         method = 'DL')

summary(res24)

plot(res24, exp = T)

# getting for 36 months - 

df36 = readxl::read_excel('F:\\GLP1_agonists\\analysis\\cox_models\\data.xlsx',
                          sheet = "36months")

glimpse(df36)

res36 = metafor::rma.uni(yi = logHR, sei = se, data = df36,
                         method = 'DL')

summary(res36)

# getting for 48 months - 

df48 = readxl::read_excel('F:\\GLP1_agonists\\analysis\\cox_models\\data.xlsx',
                          sheet = "48months")

glimpse(df48)

res48 = metafor::rma.uni(yi = logHR, sei = se, data = df48,
                         method = 'DL')

summary(res48)


