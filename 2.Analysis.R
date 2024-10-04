library(readxl)
library(ggplot2)
library(lme4)
library(patchwork)
library(reshape2)

Final_analysis_input <- 
  read.csv("Input/Final_analysis_input.csv", sep=";")

probeNames <- 
  colnames(Final_analysis_input)[startsWith(colnames(Final_analysis_input), "cg")]

Table_2 <- data.frame(
  probeName = probeNames,
  Mean_beta_value = NA,
  SD_beta_value = NA,
  Spearman_cor_estimate = NA,
  Spearman_p.value = NA
)

for (i in 1:nrow(Table_2)){
  
  spearman_input <- 
    Final_analysis_input[,c("AUC_norm",probeNames[i])]
  
  Table_2[i,c("Mean_beta_value","SD_beta_value")] <-
    c(mean(spearman_input[,probeNames[i]]), sd(spearman_input[,probeNames[i]]))
  
  spearman_result <- 
    cor.test(x = spearman_input[,probeNames[i]], y = spearman_input$AUC_norm, method = "spearman")
  
  
  Table_2[i,"Spearman_cor_estimate"] <- as.numeric(spearman_result$estimate)
  Table_2[i,"Spearman_p.value"] <- spearman_result$p.value
  
}


Table_2 <- Table_2[order(Table_2$probeName),]
Table_2$Mean_beta_value <- signif(Table_2$Mean_beta_value,2)
Table_2$SD_beta_value <- signif(Table_2$SD_beta_value,2)
Table_2$Spearman_cor_estimate <- signif(Table_2$Spearman_cor_estimate,2)

write.table(Table_2,  file = paste0("Output/Table_2.csv"), row.names = F, sep = ";")


plotList <- list()

for (probeName in Table_2$probeName){
  
  tmp <- ggplot(data = Final_analysis_input, aes_string(x = probeName, y = "AUC_norm")) + 
    geom_point() + geom_smooth(method = "lm") + labs(
      x = paste0("Beta value probe ",probeName),
      y = expression({AUC['0-24h,norm']})
    ) + theme_minimal()
  
  plotList[[probeName]] <- tmp
  
  ggsave(
    filename = paste0("Output/Figures/Exposure_vs_methylation_",probeName,".png"),
    plot = tmp,
    device = "png",
    dpi = "retina"
  ) 
  
}

plotList <- plotList[names(plotList) != "cg19046783"]

plotList[[1]] <- plotList[[1]] + ylab("")
plotList[[2]] <- plotList[[2]] + ylab("")
plotList[[3]] <- plotList[[3]] + ylab("")
plotList[[5]] <- plotList[[5]] + ylab("")
plotList[[6]] <- plotList[[6]] + ylab("")
plotList[[7]] <- plotList[[7]] + ylab("")
plotList[[8]] <- plotList[[8]] + ylab("")
plotList[[9]] <- plotList[[9]] + ylab("")

combined_plot <-  (plotList[[1]] + plotList[[2]] + plotList[[3]]) /
                  (plotList[[4]] + plotList[[5]] + plotList[[6]]) / 
                  (plotList[[7]] + plotList[[8]] + plotList[[9]])

ggsave(
  filename = paste0("Output/Figure_S2.png"),
  plot = combined_plot,
  device = "png",
  dpi = "retina"
) 

cor.test(x = Final_analysis_input$cg19046783, y = Final_analysis_input$AUC_norm, method = "pearson")

# Covariate test results (Table S1): 
summary(lm(AUC_norm ~ HCT, data = Final_analysis_input))
summary(lm(AUC_norm ~ log(ALAT), data = Final_analysis_input))
summary(lm(AUC_norm ~ CRP_above_49, data = Final_analysis_input))
summary(lm(AUC_norm ~ Recipient_CYP3A4_star22_rs35599367, data = Final_analysis_input))
summary(lm(AUC_norm ~ Donor_CYP3A4_star22_rs35599367, data = Final_analysis_input))