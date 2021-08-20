devtools::install_github("mcguinlu/robvis")

library(robvis)
library(ggplot2)
library(tidyverse)
library(readxl)

# get the dataaset

df <- read_excel(
  "F:\\GLP1_agonists\\rob\\rob_data.xlsx"
)
# plot 

summary_rob <- rob_summary(
  data = df, tool = "ROB2"
)

summary_rob

# plot traffic light

tf <- robvis::rob_traffic_light(data = df,
                                color = "cochrane",
                                tool = "ROB2")

str(tf)

ggsave(tf,
       filename = "F:\\GLP1_agonists\\rob\\rob_plot.pdf")
