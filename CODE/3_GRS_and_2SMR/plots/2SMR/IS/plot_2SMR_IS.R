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

lower_ci <- function(beta_, se_){
  #A function that takes the beta and se from summary statistics
  #and generates the lower CI.
  
  final_beta <- c()
  final_se <- c()
  
  for(index in seq(1, length(beta_))){
    
    if(as.character(beta_[index]) == "-"){
      
      final_beta <- c(final_beta, 0)
      final_se <- c(final_se, 0)
      
    } else {
      
      final_beta <- c(final_beta, beta_[index])
      final_se <- c(final_se, se_[index]) 
      
    }
    
  }
  
  lower_ci_ <- as.numeric(final_beta) - qnorm(0.975)*as.numeric(final_se)
  
}

upper_ci <- function(beta_, se_){
  #A function that takes the beta and se from summary statistics
  #and generates the upper CI.
  
  final_beta <- c()
  final_se <- c()
  
  for(index in seq(1, length(beta_))){
  
    if(as.character(beta_[index]) == "-"){
      
      final_beta <- c(final_beta, 0)
      final_se <- c(final_se, 0)
    
    } else {
    
      final_beta <- c(final_beta, beta_[index])
      final_se <- c(final_se, se_[index]) 
      
      }
    
  }
  
  upper_ci_ <- as.numeric(final_beta) + qnorm(0.975)*as.numeric(final_se)
  
  return(upper_ci_)
  
}


round_mgu_ci <- function(digit){
  
  if(is.na(digit)){
    
    return(digit)
    

  }
  
  digit_ <- round(digit, digits = 2)
  
  if(digit_ < 0){
    
    check <- abs(digit_)
    
  } else {
    
    check <- digit_
    
  }
  
  if(nchar(as.character(check)) < 4 & digit_ != 0 & digit_ != 1){
    
    #options(digits=9)
    digit_new <- as.character(paste(as.character(digit_), as.character(0), sep = ""))
    
  return(digit_new)
    
  } else {
    
    return(digit_)
    
  }
  
}

round_mgu_pval <- function(digit){
  
  if(is.na(digit)){
    
    return(digit)
    
  }
  
  if(digit < 0.01){
    
    digit_ <- formatC(digit, format = "e", digits = 2)
    
    return(digit_)
    
  } else {
    
    digit_ <- round(digit, digits = 2)
    
    return(digit_)
    
  }
  
}

clean_betas <- function(beta){
  
  if(beta == "-"){
    
    return(0)
    
  } else {
    
    return(beta)
    
  }
  
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
  nsnps <- c()
  methods <- c()
  betas <- c()
  lower_cis <- c()
  upper_cis <- c()
  ids <- c()
  pvals <- c()
  
  #First we take the exposures and methods!!
  
  exposure <- results[,1] #We can take them from here:
  method <- results[,2] #We can take them from here:
  exposure <- exposure[-1]
  method <- method[-1]
  
  exposures <- rep(exposure, 5)
  methods <- rep(method, 5)
  
  #and we get the final results without the exposures:
  
  results <- results[,-1]
  
  for(i in seq(1,25, by = 5)){
    
    j=i+4
    
    tmp=results[,seq(i,j)]
    
    #Now we are going to obtain the outcome, the beta, the lower_ci, upper_ci and pval!
    
    #First we remove the methods, cuz we do not need them anymore:
    
    tmp <- tmp[,seq(2,length(colnames(tmp)))]
    
    #First we get the outcomes, which are in the colnames
    
    outcome <- colnames(tmp)[1]
    outcome <- strsplit(outcome, "[.]")[[1]][1]
    outcome <- rep(outcome, 12)
    outcomes <- c(outcomes, outcome)
    
    #Now we remove these colnames! And get the data ready to extract the numbers:
    
    col_name_tmp <- tmp[1,]
    colnames(tmp) <- col_name_tmp
    tmp <- tmp[-1,] #removing the colnames now that we have them
    
    betas <- c(betas, as.character(unlist(sapply(as.character(tmp$b), clean_betas))))
    lower_cis <- c(lower_cis, lower_ci(tmp$b, tmp$se))
    upper_cis <- c(upper_cis, upper_ci(tmp$b, tmp$se))
    
    pvals <- c(pvals, tmp$pval)
    
    #Let's put here the number of snps:
    
    nsnps <- c(nsnps, tmp$nsnp)
    
    #And here we end the loop!
    
  }
  
  #And now we loop to get the ids
  
  for(index in seq(1, length(betas))){
    
    id <- paste(exposures[index], "-", outcomes[index], sep = "")
    
    ids <- c(ids, id)
    
  }
  
  ids <- as.data.frame(ids)
  methods <- as.data.frame(methods)
  exposures <- as.data.frame(exposures)
  outcomes <- as.data.frame(outcomes)
  nsnps <- as.data.frame(nsnps)
  betas <- as.data.frame(betas)
  lower_cis <- as.data.frame(lower_cis)
  upper_cis <- as.data.frame(upper_cis)
  pvals <- as.data.frame(pvals)
  
  final_df <- cbind(ids, methods, exposures, outcomes, betas, lower_cis, upper_cis, pvals)
  
  return(final_df)
  
}

round_mgu <- function(digit){
  
  digit_ <- round(digit, digits = 2)
  
  if(digit_ < 0){
    
    check <- abs(digit_)
    
  } else {
    
    check <- digit_
    
  }
  
  if(nchar(as.character(check)) < 4 & digit_ != 0 & digit_ != 1){
    
    #options(digits=9)
    digit_new <- as.character(paste(as.character(digit_), as.character(0), sep = ""))
    
    return(digit_new)
    
  } else {
    
    return(digit_)
    
  }
  
}

round_mgu_pval <- function(digit){
  
  if(digit < 0.01){
    
    digit_ <- formatC(digit, format = "e", digits = 2)
    
    return(digit_)
    
  } else {
    
    digit_ <- round(digit, digits = 2)
    
    return(digit_)
    
  }
  
}

##############
#Loading data#
##############

mr_results <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/supplementary_tables_per_method/GRS_&_2SMR/4_2SMR_results/2SMR_IS_results.csv")

mr_results <- as.data.frame(mr_results)

mr_results <- data_parsers(mr_results)

#################################
#Make plots for weighted section#
#################################

mr_results_og <- mr_results

mr_results_og$Effect <- paste(unlist(sapply(as.numeric(mr_results_og$betas), round_mgu)), " (", unlist(sapply(as.numeric(mr_results_og$lower_cis), round_mgu)), ";", unlist(sapply(as.numeric(mr_results_og$upper_cis), round_mgu)), ")", sep = "")

mr_results_og$alphabetical_id <- c("a1", "a2", "a3", "a4",
                                   "b1", "b2", "b3", "b4",
                                   "c1", "c2", "c3", "c4",
                                   "d1", "d2", "d3", "d4",
                                   "e1", "e2", "e3", "e4",
                                   "f1", "f2", "f3", "f4",
                                   "g1", "g2", "g3", "g4",
                                   "h1", "h2", "h3", "h4",
                                   "i1", "i2", "i3", "i4", 
                                   "j1", "j2", "j3", "j4",
                                   "k1", "k2", "k3", "k4",
                                   "l1", "l2", "l3", "l4",
                                   "m1", "m2", "m3", "m4",
                                   "n1", "n2", "n3", "n4", 
                                   "o1", "o2", "o3", "o4")                                    

mr_results$exposures <- factor(mr_results$exposures, levels = c("CIR", "CIRadjISI", "DI"))

mr_results_og$final_ids <- c("", "", "", "",
                             "", "", "FBW", "",
                             "", "", "", "",
                             "", "", "", "",
                             "", "", "FBWadjMBW", "",
                             "", "", "", "",
                             "", "", "", "",
                             "", "", "MBW", "",
                             "", "", "", "", 
                             "", "", "", "",
                             "", "", "MBWadjFBW", "",
                             "", "", "", "",
                             "", "", "", "",
                             "", "", "Father", "", 
                             "", "", "", "")

#Let's plot the ones that are OK

mr_results_og <- mr_results_og[which(!(mr_results_og$lower_cis == 0 & mr_results_og$upper_cis == 0 & mr_results_og$betas == 0)),]


plot1 <- ggplot(mr_results_og, aes(y = fct_rev(alphabetical_id), x = as.numeric(betas), shape = methods, colour = exposures)) +
  geom_point(size = 2.7, position = "identity") +  
  geom_errorbarh(aes(xmin = as.numeric(lower_cis), xmax = as.numeric(upper_cis)), height = 0.25) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  xlab("Causal effect (95% CI)") + 
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





table_base <- ggplot(mr_results_og, aes(y=alphabetical_id)) +
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


## 95% CI table
tab1 <- table_base +
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = Effect), size = 2.7,lineheight = 1,
            inherit.aes = TRUE,vjust = 1)

tab2 <- table_base +
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = unlist(sapply(as.numeric(pvals), round_mgu_pval))), size = 2.7,lineheight = 1,
            inherit.aes = TRUE, vjust = 1)

lay <-  matrix(c(1,1,1,2,2,2,2,2,2,3,3,3,4,4,4,4,4,4,4,4), nrow = 1)

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/figures/GC/IS_GC_plot.tiff", width = 3000, height=1200, res = 300)
grid.arrange(tab0, tab1, tab2, plot1, layout_matrix = lay)
dev.off()
