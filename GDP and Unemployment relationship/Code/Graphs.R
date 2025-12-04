####
#Partie 2 : Analysing the data and including it in a document
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#!!!! Make sure to change setwd() in the right depository !!!!!!

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#---------------------------------------------------------------
####

setwd("C:\\Users\\Ribak\\Documents\\ProjetPersonnel\\GDP and Unemployment relationship")  # Make sure to change to change in the right depository


rm(list = ls())

library(readxl)
library(dplyr)
library(tidyr)
#install.packages("writexl")    # Please install those packages by removing the # and run the command
library(writexl)
#install.packages("openxlsx")    # Please install those packages by removing the # and run the command
library(openxlsx)
#install.packages("ggplot2")     # Please install those packages by removing the # and run the command
library(ggplot2)



### Part 2) Analysing the data and including it in a document

# Load database, compute average for each country and compute the median

load("Data\\inflation.Rdata")

data_2013 <- subset(data_EU_inflation, as.numeric(year) >= 2013 & as.numeric(year) <= 2022)

inflation_avg_country <- data_2013 %>% 
  group_by(Country) %>%
  summarize(avg_inflation = mean(value, na.rm = TRUE))

inflation_median <- median(inflation_avg_country$avg_inflation, na.rm = TRUE)

countries_above_median_inflation <- inflation_avg_country %>%
  filter(avg_inflation > inflation_median) %>%
  select(Country)


# Create table and a bar chart

bar_chart <- ggplot(inflation_avg_country, aes(x = reorder(Country, avg_inflation), y = avg_inflation)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = inflation_median, linetype = "dashed", color = "red", size = 1) +
  coord_flip() +
  labs(
    title = "Average Inflation Rate (2013–2022) by Country",
    subtitle = "EU27 Median Inflation Shown as a Red Dashed Line",
    x = "Country",
    y = "Average Inflation Rate (%)"
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 12)
  )

ggsave("Graph\\European_inflation_chart_bar.png", plot = bar_chart, width = 10, height = 10)


#  list all the countries with an average inflation above the median as a table in your document

countries_above_median_inflation <- inflation_avg_country %>%
  filter(avg_inflation > inflation_median) %>%
  select(Country)

list_countries <- countries_above_median_inflation %>%
  group_by(Country) %>%
  summarize(
    "Inflation average"  = "Above")
    

print(list_countries)  

write.xlsx(list_countries, "Data\\list_country.xlsx")


# List all countries and there average inflation. Compute the median inflation 

list_avg_inflation <- inflation_avg_country %>%
  mutate(median = inflation_median) %>%
  select(Country, inflation = avg_inflation, median)

write.xlsx(list_avg_inflation, "Data\\List of countries with average and median inflation.xlsx")
  
  #################################################################################################################################


## Load unemployment data, compute average for each country and his median

load("Data\\unemployment.RData")
summary(data_2013)

# Change year in the right format (numeric)

data_EU_unemployment$year <- as.numeric(data_EU_unemployment$year)

# Filter the data the have above and below label for each countries concerned by the filter

unemployment_2013_2022 <- data_EU_unemployment %>%
  filter(year >= as.numeric(2013) & year <= as.numeric(2022)) %>%
  select(Country,year,value)


unemployment_2013_2022_above_median <- unemployment_2013_2022 %>%
  filter(Country %in% countries_above_median_inflation$Country) %>%
  select(Country,year,value)



# Create a plot of the unemployment rate development


line_chat <- ggplot(unemployment_2013_2022_above_median, aes(x = year, y = value, color = Country , group = Country)) +
  geom_line(size = 1) +
  labs(
    title = "Unemployment rate between 2013-2022",
    x = "Year",
    y = "Unemployment rate in %",
    color = "Country"
     ) +
  theme_minimal()+
  theme(
    axis.title = element_text(size = 14),
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 22),
    axis.line = element_line(color = "black"),
    panel.background = element_rect( fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, size = 28, margin = margin(t = 20,b = 25)),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90")
  )

ggsave("Graph\\Unemployment_rate_plot_above.png", plot = line_chat, width = 15, height = 10, bg = "white")


###########  Draw aboxplot diagram of each country in a single plot. ###################################################


boxplot_chart <- ggplot(unemployment_2013_2022_above_median, aes(x = Country, y = value)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  labs(
    title = "Unemployment Rate Distribution by Country (2013-2022)",
    x = "Country",
    y = "Unemployment Rate (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate country names
    axis.title = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 22, margin = margin(t = 20, b = 25)),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90")
  )

ggsave("Graph\\unemployment_boxplot.png", plot = boxplot_chart, width = 15, height = 10, bg = "white")


########################################################################################################
########################################################################################################


data_EU_inflation_and_unemployment <- merge(data_EU_inflation, data_EU_unemployment, by = c("Country", "year"), all = TRUE, na.rm = TRUE)

data_EU_inflation_and_unemployment$year <- as.numeric(data_EU_inflation_and_unemployment$year)


data_EU_inflation_and_unemployment <- data_EU_inflation_and_unemployment %>%
  filter(year >= as.numeric(2013) & year <= as.numeric(2022))


##   Generate a variable with the first three uppercase letters of each country’s name

data_EU_inflation_and_unemployment$country_initials <- toupper(substr(data_EU_inflation_and_unemployment$Country,1,3))

unique(data_EU_inflation_and_unemployment$country_initials)

##   Create a scatter plot showing the relationship between average unemployment
# Create data such as average inflation and average unemployment

average_data_inflation_unemployment <- data_EU_inflation_and_unemployment %>%
  group_by(Country) %>%
  summarize(
    avg_inflation = mean(value.x, na.rm = TRUE),
    avg_unemployment = mean(value.y, na.rm = TRUE)
  )

average_data_inflation_unemployment$Country_initials <- toupper(substr(average_data_inflation_unemployment$Country, 1 , 3))

# Create the scatter plot

scatter_plot <- ggplot(average_data_inflation_unemployment, aes(x = avg_unemployment, y = avg_inflation, group = Country))+
  geom_point(color = "blue", size = 3) +
  geom_text(aes (label = Country_initials), vjust = -1, hjust = 0.5, size = 4) +
  labs(
    title = "Relationship Between Average Unemployment and Average Inflation (2013-2022)",
    x = "Average Unemployment Rate (%)",
    y = "Average Inflation Rate (%)"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 22, margin = margin(t = 20, b = 25)),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90")
  )

ggsave("Graph\\relation_unemployment_inflation_reversed.png", width = 15, height = 10, bg = "white")  
  
  
  #######################################################################################################################
#########################################################################################################################

# Load your GDP data set and merge it with your unemployment data set.

load("Data\\GDP.RData")

#Normalise data for the merge
data_EU_GDP$Country <- trimws(tolower(data_EU_GDP$Country))
data_EU_unemployment$Country <- trimws(tolower(data_EU_unemployment$Country))
data_EU_GDP$year <- as.numeric(data_EU_GDP$year)
data_EU_unemployment$year <- as.numeric(data_EU_unemployment$year)

data_GDP_unemployment_merged <- merge(data_EU_GDP, data_EU_unemployment, by = c("Country", "year")) 

# Create a new variable of log GDP 

data_GDP_unemployment_merged$lgdp <- log(data_GDP_unemployment_merged$value.x)

table_statistics <- data_GDP_unemployment_merged %>%
  group_by(year) %>%
  summarize(
    avg_unemployment = mean(value.y, na.rm = TRUE),
    avg_gdp = mean(value.x, na.rm = TRUE),
    min_unemployment = min(value.y, na.rm = TRUE),
    max_unemployment = max(value.y, na.rm = TRUE),
    median_unemployment = median(value.y, na.rm = TRUE),
    min_gdp = min(value.x, na.rm = TRUE),
    max_gdp = max(value.x, na.rm = TRUE),
    median_gdp = median(value.x, na.rm = TRUE),
    
  )

print(table_statistics)  

write.xlsx(table_statistics, "Data\\Tables_statistique.xlsx") 


# Create a line plot showing Belgium’s logged GDP on one y-axis over time.
## filter data to choose only Belgium 

data_belgium <- data_GDP_unemployment_merged %>%
  filter(Country == "belgium")


# Plot the graph

pdf("Graph\\Log gdp and unemployment graph.pdf")

# Adjust the margin to fit everything in the graph
par(mar = c(13,4,4,5)+ 0.1)

# Plot the first line which is log GDP
plot(data_belgium$year,
     data_belgium$lgdp,
     col = "blue",
     type = "l",
     ylab = "GDP",
     xlab= "year",
     lwd = 1.5)

# Fit the other graph into the same one

par(new = TRUE)

# Plot the second line which is unemployment

plot(data_belgium$year,
     data_belgium$value.y,
     col = "red",
     type = "l",
     ylab = "",
     xlab= "",
     lwd = 1.5,
     axes = FALSE)

# add text for the second axis

axis(4 , col = "black")

mtext("Unemployment",
      side = 4,
      line = 3)

# add a legend
legend("bottom",
       legend = c("Log GDP", "Unemployment"),
       col = c("blue", "red"),
       lty = 1,
       lwd = 1.5,
       bty = "n",
       cex = 1.1,
       xpd = TRUE,
       horiz = TRUE,
       inset = c(0, -0.20),
       xjust = 0.5)

dev.off()





