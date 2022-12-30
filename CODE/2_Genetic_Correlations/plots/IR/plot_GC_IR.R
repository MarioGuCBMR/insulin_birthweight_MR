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
  
  exposure <- results$row
  exposure <- exposure[-1]
  exposures <- rep(exposure, 5)
  
  #and we get the final results without the exposures:
  
  results <- results[,-1]
  
  for(i in seq(1,20, by = 4)){
    
    j=i+3
    
    tmp=results[,seq(i,j)]
    
    #Now we are going to obtain the outcome, the beta, the lower_ci, upper_ci and pval!
    
    #First we get the outcomes, which are in the colnames
    
    outcome <- colnames(tmp)[1]
    outcome <- rep(outcome, 3)
    outcomes <- c(outcomes, outcome)
    
    #Now we remove these colnames! And get the data ready to extract the numbers:
    
    col_name_tmp <- tmp[1,]
    colnames(tmp) <- col_name_tmp
    tmp <- tmp[-1,] #removing the colnames now that we have them
    
    betas <- c(betas, tmp$`Genetic Correlation`)
    lower_cis <- c(lower_cis, tmp$`95% LCI`)
    upper_cis <- c(upper_cis, tmp$`95% UCI`)
    
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

round_mgu <- function(digit){
  
  digit_ <- round(digit, digits = 2)
  
  if(digit_ < 0){
    
    check <- abs(digit_)
    
  } else {
    
    check <- digit_
    
  }
  
  if(nchar(as.character(check)) < 4){
    
    #options(digits=9)
    digit_new <- as.character(paste(as.character(digit_), as.character(0), sep = ""))
    
    return(digit_new)
    
  } else {
    
    return(digit_)
    
  }
  
}

##############
#Loading data#
##############

gc_results <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/supplementary_tables_per_method/GC/results_IR_GC.csv")

gc_results <- as.data.frame(gc_results)

gc_results <- data_parsers(gc_results)

#################################
#Make plots for weighted section#
#################################

gc_results_og <- gc_results

gc_results_og$CI <- paste("(", unlist(sapply(as.numeric(gc_results_og$lower_cis), round_mgu)), ";", unlist(sapply(as.numeric(gc_results_og$upper_cis), round_mgu)), ")", sep = "")

gc_results_og$alphabetical_id <- c("a1", "a2", "a3", "a4", "a5",
                                         "b1", "b2", "b3", "b4", "b5",
                                         "c1", "c2", "c3", "c4", "c5")

plot1 <- ggplot(gc_results_og, aes(y = fct_rev(alphabetical_id), x = as.numeric(betas), shape = exposures)) +
  geom_point(size = 2.7, position = "identity") +  
  geom_errorbarh(aes(xmin = as.numeric(lower_cis), xmax = as.numeric(upper_cis)), height = 0.25) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  xlab("GC (95% CI)") + 
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



gc_results_og$final_ids <- c("", "FBW", "",
                                   "", "FBWadjMBW", "",
                                   "", "MBW", "",
                                   "", "MBWadjFBW", "",
                                   "", "Father", "")

table_base <- ggplot(gc_results_og, aes(y=alphabetical_id)) +
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
            inherit.aes = TRUE, vjust = 1)


tab1 <- table_base + 
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = unlist(sapply(as.numeric(betas), round_mgu))), size =2.7,lineheight = 1,
            inherit.aes = TRUE, vjust = 1)


## 95% CI table
tab2 <- table_base +
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = CI), size = 2.7,lineheight = 1,
            inherit.aes = TRUE,vjust = 1)

tab3 <- table_base +
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = pvals), size = 2.7,lineheight = 1,
            inherit.aes = TRUE, vjust = 1)

lay <-  matrix(c(1,1,1,2,2,3,3,3,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5), nrow = 1)

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/figures/GC/IR_GC_plot.tiff", width = 3000, height=1200, res = 300)
grid.arrange(tab0, tab1, tab2, tab3,plot1, layout_matrix = lay)
dev.off()
