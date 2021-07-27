# forest plot for myocardial infarction
# run the parametric model and then combine to plot.
# INCOMPLETE...

library(easypackages)
libraries(c("metaRMST","patchwork","tidyverse",
            "ckbplotr"))


sustain <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/mi_sustain6.csv')

sustain$trialID <- 1

harmony <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/mi_harmony.csv')

harmony$trialID <- 2

leader <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/mi_leader.csv')

leader$trialID <- 3

rewind <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/mi_rewind.csv')

rewind$trialID <- 4



df <- rbind(rewind, leader, sustain, 
            harmony)


df$triaID <- df$trialID
df$Time <- df$time
df$Event <- df$status
df$Arm <- df$treat

df2 <- df %>% select(trialID, Time, Event, Arm)

# fitting the parametric model 


res <- metaRMSTD(df2,time_horizons = c(24,48),
                 MA_method ="uni_flex")

df_res <- tbl_df(res$result)

# obtain data for the plots.


dfm <- res$rmstd_est %>% tbl_df()

dfse <- res$se_rmstd_est %>% tbl_df()


table_24 <- tibble(
  study = c("Sustain6", 
            "Harmony Outcomes","LEADER","REWIND",NA,
            "Pooled Estimate"),
  estimate = c(dfm$RMSTD_est_at_24,NA,0.0804),
  stderr = c(dfse$se_RMSTD_est_at_24,NA, 0.0318),
  variable = c("1","2","3","4",NA,"Pooled"),
  color = c(rep("black",4),NA,"red")
)


table_48 <- tibble(
  study = c("Sustain6", 
            "Harmony Outcomes","LEADER","REWIND",NA,
            "Pooled Estimate"),
  estimate = c(dfm$RMSTD_est_at_48,NA,0.416),
  stderr = c(dfse$se_RMSTD_est_at_48,NA, 0.223),
  variable = c("1","2","3","4",NA,"Pooled"),
  color = c(rep("black",4),NA,"red")
)


rowlabels = tibble(heading = c(rep("Myocardial Infarction",6)),
                   subheading = c(rep("Trials Included",4),
                                  NA, NA),
                   variable = c("1","2","3","4",NA, "Pooled"),
                   label = c("Sustain6", 
                             "Harmony Outcomes","LEADER","REWIND",
                             NA, "Pooled Estimate"))





f <- make_forest_plot(panels = list(table_24, table_48),
                      
                      row.labels = rowlabels,
                      col.key = "variable",
                      rows = c("Myocardial Infarction"),
                      row.labels.levels = c("heading","subheading",
                                            "label"),
                      panel.names = c("dRMST at 24 months", "dRMST at 48 months"
                      ),
                      exponentiate = F,
                      col.right.heading = "RMST Diff.(95%CI)",
                      colour = "color",
                      nullval = 0,
                      xlim = c(-1,1.5),
                      xlab = "dRMST(95% CI)   favours GLP1 agonist")


f

# save this plot now.

ggsave(f$plot,
       filename = "F:\\GLP1_agonists\\analysis\\results\\mi_forest.pdf",
       height = 5,
       width = 8,
       units = "in",
       device = "pdf")


# NEED TO MAKE SOME ADJUSTMENTS AND CREATE AN ARROW ---

