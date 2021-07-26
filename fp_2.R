# create forest plot for 3 pt MACE.
# get the results from the model fits.
# combine all 4 time points into 1 panel.




test2 <- tibble(
  study = c("Amplitude-O","ELIXA","Sustain6", 
            "Harmony Outcomes",'EXCEL',"LEADER","REWIND","PIONEER6",NA,
            "Pooled Estimate"),
  estimate = c(dfm$RMSTD_est_at_48,NA,0.627),
  stderr = c(dfse$se_RMSTD_est_at_48,NA, 0.182),
  variable = c("1","2","3","4","5","6","7","8",NA,"Pooled"),
  color = c(rep("black",8),NA,"red")
)


rowlabels2 = tibble(heading = c(rep("3 point MACE",10)),
                   subheading = c(rep("Trials Included",8),
                                  NA, NA),
                   variable = c("1","2","3","4","5","6","7","8",NA, "Pooled"),
                   label = c("Amplitude-O","ELIXA","Sustain6", 
                             "Harmony Outcomes",'EXCEL',"LEADER","REWIND","PIONEER6",
                             NA, "Pooled Estimate"))


library(ckbplotr)

f2 <- make_forest_plot(panels = list(test2),
                      
                      row.labels = rowlabels2,
                      col.key = "variable",
                      rows = c("3 point MACE"),
                      row.labels.levels = c("heading","subheading",
                                            "label"),
                      panel.names = "dRMST at 48 months",
                      exponentiate = F,
                      col.right.heading = "RMST Diff.(95%CI)",
                      colour = "color",
                      nullval = 0)



f2




ggsave(plot = f2$plot,
  
  filename = "F:\\GLP1_agonists\\analysis\\results\\3ptmace_fp.pdf",
       device = "pdf",
       height = 5,
       width = 5,
       units = "in")


# making plot panels at 12, 24, 36, 48 months
# combine plot to make 1 big figure


f3 <- make_forest_plot(panels = list(test2, test2),
                       
                       row.labels = rowlabels2,
                       col.key = "variable",
                       rows = c("3 point MACE"),
                       row.labels.levels = c("heading","subheading",
                                             "label"),
                       panel.names = c("dRMST at 48 months", "dRMST at 24 months"),
                       exponentiate = F,
                       col.right.heading = "RMST Diff.(95%CI)",
                       colour = "color",
                       nullval = 0)

