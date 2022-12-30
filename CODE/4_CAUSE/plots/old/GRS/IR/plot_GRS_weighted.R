##############
#INTRODUCTION#
##############

#This code will take all the GRS results and make forest plots automatically.
#This is going to take the weighted GRS only.
#We will do the unweighted later.

###########
#Libraries#
###########

library(tidyverse)
library(ggplot2)
library(data.table)
library(gridExtra)

###################
#Loading functions#
###################

beta_parser <- function(effect_size){
  
  tmp <- as.character(strsplit(effect_size, " ")[[1]][1])
  
  return(tmp)
  
}

lower_ci_parser <- function(effect_size){
  
  tmp <- as.character(strsplit(effect_size, " ")[[1]][2])
  tmp <- as.character(strsplit(tmp, "[()]")[[1]][2])
  
  return(tmp)
  
}

upper_ci_parser <- function(effect_size){
  
  tmp <- as.character(strsplit(effect_size, " ")[[1]][4])
  tmp <- as.character(strsplit(tmp, "[)]")[[1]][1])
  
  return(tmp)
  
}

data_parsers <- function(results){
  #This function is going to return a whole dataframe with the following vectors:
  #exposure
  #outcome
  #beta
  #lower_ci
  #upper_ci
  #id (exposure-outcome)
  #pvals
  
  #First we make empty vectors to save all:
  
  outcomes <- c()
  exposures <- c()
  betas <- c()
  lower_cis <- c()
  upper_cis <- c()
  ids <- c()
  pvals <- c()
  
  #First we take the exposures!!
  
  exposure <- results[,1]
  exposure <- exposure[-1]
  exposures <- rep(exposure, 5)
  
  #and we get the final results without the exposures:
  
  results <- results[,-1]
  
  for(i in seq(1,9, by = 2)){
    
    j=i+1
    
    tmp=results[,seq(i,j)]
    
    #Now we are going to obtain the outcome, the beta, the lower_ci, upper_ci and pval!
    
    #First we get the outcomes, which are in the colnames
    
    outcome <- colnames(tmp)[1]
    outcome <- rep(outcome, 6)
    outcomes <- c(outcomes, outcome)
    
    #Now we remove these colnames! And get the data ready to extract the numbers:
    
    col_name_tmp <- tmp[1,]
    colnames(tmp) <- col_name_tmp
    tmp <- tmp[-1,] #removing the colnames now that we have them
    
    beta <- unlist(sapply(tmp$`GRS effect size (95% CI)`, beta_parser))
    lower_ci <- unlist(sapply(tmp$`GRS effect size (95% CI)`, lower_ci_parser))
    upper_ci <- unlist(sapply(tmp$`GRS effect size (95% CI)`, upper_ci_parser))
    
    betas <- c(betas, beta)
    lower_cis <- c(lower_cis, lower_ci)
    upper_cis <- c(upper_cis, upper_ci)
    
    pvals <- c(pvals, tmp$`P-value`)
    
    #And here we end the loop!
    
  }
  
  #And now we loop to get the ids
  
  for(index in seq(1, length(betas))){
    
    id <- paste(exposures[index], "-", outcomes[index], sep = "")
    
    ids <- c(ids, id)
    
  }
  
  ids <- as.data.frame(ids)
  exposures <- as.data.frame(exposures)
  outcomes <- as.data.frame(outcomes)
  betas <- as.data.frame(betas)
  lower_cis <- as.data.frame(lower_cis)
  upper_cis <- as.data.frame(upper_cis)
  pvals <- as.data.frame(pvals)
  
  final_df <- cbind(ids, exposures, outcomes, betas, lower_cis, upper_cis, pvals)
  
  return(final_df)
  
}

##############
#Loading data#
##############

weighted_results <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/supplementary_tables_per_method/GRS_&_2SMR/3_results_4_GRS/results_GRS_weighted.csv")
unweighted_results <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/supplementary_tables_per_method/GRS_&_2SMR/3_results_4_GRS/results_GRS_unweighted.csv")

weighted_results <- as.data.frame(weighted_results)
unweighted_results <- as.data.frame(unweighted_results)

weighted_results <- data_parsers(weighted_results)
unweighted_results <- data_parsers(unweighted_results)

#################################
#Make plots for weighted section#
#################################

weighted_results_og <- weighted_results

weighted_results_og$CI <- paste("(", weighted_results_og$lower_cis, ";", weighted_results_og$upper_cis, ")", sep = "")

weighted_results_og$alphabetical_id <- c("a1", "a2", "a3", "a4", "a5",
                                         "b1", "b2", "b3", "b4", "b5",
                                         "c1", "c2", "c3", "c4", "c5",
                                         "d1", "d2", "d3", "d4", "d5",
                                         "e1", "e2", "e3", "e4", "e5", 
                                         "f1", "f2", "f3", "f4", "f5")

plot1 <- ggplot(weighted_results_og, aes(y = fct_rev(alphabetical_id), x = as.numeric(betas), col=outcomes, shape = exposures)) +
  geom_point(size = 2.7, position = "identity") +  
  geom_errorbarh(aes(xmin = as.numeric(lower_cis), xmax = as.numeric(upper_cis)), height = 0.25) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  xlab("GRS (95% CI)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_blank(),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.ticks.y = element_blank(),
        legend.position = "none")



weighted_results_og$final_ids <- c("", "FBW", "", "", "", "",
                                   "", "FBWadjMBW", "", "", "", "",
                                   "", "MBW", "", "", "", "",
                                   "", "MBWadjFBW", "", "", "", "",
                                   "", "Father", "", "", "", "")

table_base <- ggplot(weighted_results_og, aes(y=alphabetical_id)) +
  ylab(NULL) + xlab("  ") + 
  theme(plot.title = element_text(hjust = 0.5, size=12), 
        axis.text.x = element_text(color="white", hjust = -3, size = 25), ## This is used to help with alignment
        axis.line = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),
        axis.title.y = element_blank(), 
        legend.position = "none",
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.background = element_blank())

tab0 <- table_base + 
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = final_ids), size = 3, lineheight = .2,
            inherit.aes = TRUE, position = position_nudge(y = -0.5), vjust = 1)


tab1 <- table_base + 
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = betas), size =2.7,lineheight = 1,
            inherit.aes = TRUE, vjust = 1)


## 95% CI table
tab2 <- table_base +
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = CI), size = 2.7,lineheight = 1,
            inherit.aes = TRUE,vjust = 1)

tab3 <- table_base +
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = pvals), size = 2.7,lineheight = 1,
            inherit.aes = TRUE, vjust = 1)


#################################################
#NOW LET'S MAKE THE DATA FOR THE UNWEIGHTED ONES#
#################################################

unweighted_results_og <- unweighted_results

unweighted_results_og$CI <- paste("(", unweighted_results_og$lower_cis, ";", unweighted_results_og$upper_cis, ")", sep = "")

unweighted_results_og$alphabetical_id <- c("a1", "a2", "a3", "a4", "a5",
                                         "b1", "b2", "b3", "b4", "b5",
                                         "c1", "c2", "c3", "c4", "c5",
                                         "d1", "d2", "d3", "d4", "d5",
                                         "e1", "e2", "e3", "e4", "e5", 
                                         "f1", "f2", "f3", "f4", "f5")

plot2 <- ggplot(unweighted_results_og, aes(y = fct_rev(alphabetical_id), x = as.numeric(betas), col=outcomes, shape = exposures)) +
  geom_point(size = 2.7, position = "identity") +  
  geom_errorbarh(aes(xmin = as.numeric(lower_cis), xmax = as.numeric(upper_cis)), height = 0.25) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  xlab("GRS (95% CI)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_blank(),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.ticks.y = element_blank(),
        legend.position = "right")



unweighted_results_og$final_ids <- c("", "FBW", "", "", "", "",
                                     "", "FBWadjMBW", "", "", "", "",
                                     "", "MBW", "", "", "", "",
                                     "", "MBWadjFBW", "", "", "", "",
                                     "", "Father", "", "", "", "")

table_base.2 <- ggplot(unweighted_results_og, aes(y=alphabetical_id)) +
  ylab(NULL) + xlab("  ") + 
  theme(plot.title = element_text(hjust = 0.5, size=12), 
        axis.text.x = element_text(color="white", hjust = -3, size = 25), ## This is used to help with alignment
        axis.line = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),
        axis.title.y = element_blank(), 
        legend.position = "none",
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.background = element_blank())

tab0.2 <- table_base.2 + 
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = final_ids), size = 3, lineheight = .2,
            inherit.aes = TRUE, position = position_nudge(y = -0.5), vjust = 1)


tab1.2 <- table_base.2 + 
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = betas), size =2.7,lineheight = 1,
            inherit.aes = TRUE, vjust = 1)


## 95% CI table
tab2.2 <- table_base.2 +
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = CI), size = 2.7,lineheight = 1,
            inherit.aes = TRUE,vjust = 1)

tab3.2 <- table_base.2 +
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = pvals), size = 2.7,lineheight = 1,
            inherit.aes = TRUE, vjust = 1)



lay <-  matrix(c(1,1,1,2,2,3,3,3,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,7,7,8,8,8,8,8,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10), nrow = 1)
plot_end <- grid.arrange(tab0, tab1, tab2, tab3,plot1, tab0.2, tab1.2, tab2.2, tab3.2,plot2, layout_matrix = lay)

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/figures/GRS/GRS_plot.tiff", width = 5000, height=3000, res = 300)
grid.arrange(tab0, tab1, tab2, tab3,plot1, tab0.2, tab1.2, tab2.2, tab3.2,plot2, layout_matrix = lay)
dev.off()
