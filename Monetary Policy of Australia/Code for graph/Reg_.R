setwd("C:\Users\\Ribak\\Documents\\ProjetPersonnel\\Monetary Policy of Australia\\Data")


rm(list = ls())

data <- read_excel("Data.xlsx")

model <- lm(interest_rate ~ Outputgap + Inflation + Lag_IR, data = data)

summary(model)

model_before_2008 <- subset(data, data$Year <= 2008)
model_after_2008 <- subset(data, data$Year > 2008)


model_avant_2008 <- lm(interest_rate ~ Outputgap + Inflation + Lag_IR,
                       model_before_2008 = model_before_2008)

model_apres_2008 <- lm(interest_rate ~ Outputgap + Inflation + Lag_IR,
                       model_before_2008 = model_after_2008)
