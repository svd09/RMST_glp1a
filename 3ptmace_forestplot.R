# create forest plot for 3 pt MACE.
# get the results from the model fits.
# combine all 4 time points into 1 panel.



library(easypackages)

libraries(c("tidyverse","IPDfromKM",'survival',
            "flexsurv","broom","rstpm2","survminer",
            "metaRMST", "ggthemes","ckbplotr"))



table_12 <- read_csv("F:/GLP1_agonists/data/forest_plot_data/3ptmace/table_12.csv")

table_24 <- read_csv("F:/GLP1_agonists/data/forest_plot_data/3ptmace/table_24.csv")

table_36 <- read_csv("F:/GLP1_agonists/data/forest_plot_data/3ptmace/table_36.csv")
                     
table_48 <- read_csv("F:/GLP1_agonists/data/forest_plot_data/3ptmace/table_48.csv")  


rowlabels <- read_csv("F:/GLP1_agonists/data/forest_plot_data/3ptmace/3ptmace_rowlabels2.csv")



f <- make_forest_plot(panels = list(table_12, table_24),
                       
                       row.labels = rowlabels,
                       col.key = "variable",
                       rows = c("3 point MACE"),
                       row.labels.levels = c("heading","subheading",
                                             "label"),
                       panel.names = c("dRMST at 12 months", "dRMST at 24 months"
                                       ),
                       exponentiate = F,
                       col.right.heading = "RMST Diff.(95%CI)",
                       colour = "color",
                       nullval = 0)


f


f2 <- make_forest_plot(panels = list(table_36, table_48),
                       
                       row.labels = rowlabels,
                       col.key = "variable",
                       rows = c("3 point MACE"),
                       row.labels.levels = c("heading","subheading",
                                             "label"),
                       panel.names = c("dRMST at 36 months", "dRMST at 48 months"
                       ),
                       exponentiate = F,
                       col.right.heading = "RMST Diff.(95%CI)",
                       colour = "color",
                       nullval = 0)

f2