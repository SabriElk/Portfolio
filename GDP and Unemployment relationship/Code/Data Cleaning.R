####
#Partie 1 : Create a data set with EU countries
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#!!!! Make sure to change setwd() in the right depository !!!!!!

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#---------------------------------------------------------------
#### This code has been written by Sabri El kaddouri

setwd("C:\\Users\\Ribak\\Documents\\ProjetPersonnel\\GDP and Unemployment relationship") # Make sure to change to change in the right depository

#install.packages("dplyr")
#install.packages("readxl")
#install.packages("tidyr")
library(readxl)
library(dplyr)
library(tidyr)

rm(list = ls())

#  Setup the data

data_country <- read_excel("RawData\\Country_Codes_and_Names.xlsx")
drop(data_country$AREA)

#  Keeping only the country part of the European Union, excluding United Kingdom

Eu_members <- subset(data_country, data_country$AREA == "European Union (EU)")
Eu_members <- subset(Eu_members, Eu_members$`COUNTRY NAME` != "United Kingdom")

#  Rename countries the right. In that case rename Czech Republic into Czechia
# And Germany (including former GDR)

Eu_members$`COUNTRY NAME`[Eu_members$`COUNTRY NAME` == "Czech Republic"] <-  "Czechia"
Eu_members$`COUNTRY NAME`[Eu_members$`COUNTRY NAME` == "Germany (including former GDR from 1991)"] <-  "Germany"


#  Save Eumember data set in a file

save(Eu_members, file = "Data\\Eu_members.Rdata")



#############################################################################################

#                  Partie 1.2 : Import and reshape panel data                               #

#############################################################################################

?read.csv

# a)  Import the data set.

data_inflation <- read.csv("RawData\\API_FP.CPI.TOTL.ZG_DS2_en_csv_v2_77.csv",
                           header = FALSE)
data_GDP <- read_excel("RawData\\nama_10_gdp$defaultview_page_spreadsheet.xlsx",
                       sheet = 3,
                       range = "A9:K54")

data_unemployment <- read.csv("RawData\\API_SL.UEM.TOTL.ZS_DS2_en_csv_v2_41.csv",
                          header = FALSE)

# b) Rename the variables so that each indicator is named in the following way: indicator_year

start_year <- 1960
year <- seq(start_year, start_year + 69-5)
colnames(data_inflation)[5:69] <- paste("inflation_",year)


colnames(data_unemployment)[5:69] <- paste("unemployment_", year)

start_year_2014 <- 2014
year_2014 <- seq(start_year_2014, start_year_2014 + 11-2)
colnames(data_GDP)[2:11] <- paste("GDP_", year_2014)



#  c)  Make sure all numeric variables are also saved as numeric variables. This might imply some cleaning

data_GDP <- data_GDP %>% mutate(across(starts_with("GDP_"), as.numeric))
data_inflation <- data_inflation %>% mutate(across(starts_with("inflation_"), as.numeric))
data_unemployment <- data_unemployment %>% mutate(across(starts_with("unemployment"), as.numeric))


data_inflation$`V1`[data_inflation$`V1` == "Slovak Republic"] <- "Slovakia"
data_unemployment$`V1`[data_unemployment$`V1` == "Slovak Republic"] <- "Slovakia"

data_country[data_country == ":"] <- NA 
data_GDP[data_GDP == ":"] <- NA 
data_inflation[data_inflation == ":"] <- NA 
data_unemployment[data_unemployment == ":"] <- NA 


#  d)  Merge the "EU_members" data to each of the three data sets
# Rename variable of interest to Country

colnames(data_GDP)[colnames(data_GDP) == "TIME"] <- "Country"
colnames(data_inflation)[colnames(data_inflation) == "V1"] <- "Country"
colnames(data_unemployment)[colnames(data_unemployment) == "V1"] <- "Country"
colnames(Eu_members)[colnames(Eu_members) == "COUNTRY NAME"] <- "Country"

# Merge each data set with EU_members

data_EU_GDP <- merge(Eu_members, data_GDP, by = "Country")
data_EU_unemployment <- merge(Eu_members, data_unemployment, by = "Country")
data_EU_inflation  <- merge(Eu_members, data_inflation, by = "Country")

# e)  Reshape the data set to be in the long format
# Drop all variable useless

data_EU_GDP$AREA <- NULL
data_EU_GDP$CODE <- NULL
data_EU_GDP$AREA <- NULL

data_EU_inflation$AREA <- NULL
data_EU_inflation$CODE <- NULL
data_EU_inflation$V2 <- NULL
data_EU_inflation$V3 <- NULL
data_EU_inflation$V4 <- NULL

data_EU_unemployment$AREA <- NULL
data_EU_unemployment$CODE <- NULL
data_EU_unemployment$V2 <- NULL
data_EU_unemployment$V3 <- NULL
data_EU_unemployment$V4 <- NULL

# Reshape in long format 

data_EU_GDP <- data_EU_GDP %>%
  pivot_longer(
    cols = -Country,
    names_to = c("indicator","year"),
    names_sep = "_"
  )

data_EU_inflation <- data_EU_inflation %>% 
  pivot_longer(
    cols = -Country,
    names_to = c("indicator","year"),
    names_sep = "_"
  )

data_EU_unemployment <- data_EU_unemployment %>% 
  pivot_longer(
    cols = -Country,
    names_to = c("indicator","year"),
    names_sep = "_"
  )


#  Save the data set named accordingly to the respective indicator


data_EU_GDP$indicator <- NULL
data_EU_inflation$indicator <- NULL
data_EU_unemployment$indicator <- NULL




save(data_EU_GDP, file = "Data\\GDP.Rdata")
save(data_EU_inflation, file = "Data\\inflation.Rdata")
save(data_EU_unemployment, file = "Data\\unemployment.Rdata")



