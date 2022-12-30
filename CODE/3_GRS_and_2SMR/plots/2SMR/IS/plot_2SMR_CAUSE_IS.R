##############
#INTRODUCTION#
##############

#This code will take all the 2SMR and CAUSE results and make forest plots automatically.

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
  
  if(is.na(digit)){
    
    digit_ <- ""
    
    return(digit_)
    
  }
  
  if(digit < 0.01){
    
    digit_ <- formatC(digit, format = "e", digits = 2)
    
    return(digit_)
    
  } else {
    
    digit_ <- round(digit, digits = 2)
    
    return(digit_)
    
  }
  
}

round_mgu_cause <- function(digit){
  
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

beta_parser <- function(effect_size){
  
  tmp <- as.character(strsplit(effect_size, " ")[[1]][1])
  
  return(tmp)
  
}

lower_ci_parser_cause <- function(effect_size){
  
  tmp <- as.character(strsplit(effect_size, " ")[[1]][2])
  tmp <- as.character(strsplit(tmp, "[()]")[[1]][2])
  tmp <- as.character(strsplit(tmp, ",")[[1]][1])
  
  
  return(tmp)
  
}

upper_ci_parser_cause <- function(effect_size){
  
  tmp <- as.character(strsplit(effect_size, " ")[[1]][3])
  tmp <- as.character(strsplit(tmp, "[)]")[[1]][1])
  
  return(tmp)
  
}

model_parser <- function(tmp_){
  #This gets the cause results for a exposure-outcome in their raw version
  #and cleans them.
  #So that we select the causal/shared effect in their correct format.
  
  #And now we are going to select the model that works for them:
  
  index_ <- ifelse(tmp_$p[1] < 0.05, 2, 1) #select the second row if significant
  
  if(index_ == 2){ #if it chooses the causal model...
    
    tmp_end <- tmp_[index_,]
    
    tmp_end_clean <- tmp_end %>%
      select(model, gamma, q, p)
    
    colnames(tmp_end_clean) <- c("model", "beta", "q", "p")
    
  } else {
    
    tmp_end <- tmp_[index_,]
    
    tmp_end_clean <- tmp_end %>%
      select(model, eta, q, p)
    
    colnames(tmp_end_clean) <- c("model", "beta", "q", "p")
    
  }
  
  #Now we just need to clean the beta, which is such a pain to do.
  #We will do it with our cleaning functions :)
  
  beta_effect_clean <- beta_parser(tmp_end_clean$beta)
  beta_lower_ci_clean <- lower_ci_parser_cause(tmp_end_clean$beta)
  beta_upper_ci_clean <- upper_ci_parser_cause(tmp_end_clean$beta)
  
  beta_clean_end <- paste(as.character(round_mgu_cause(as.numeric(beta_effect_clean))),
                          " (",
                          as.character(round_mgu_cause(as.numeric(beta_lower_ci_clean))),
                          ", ",
                          as.character(round_mgu_cause(as.numeric(beta_upper_ci_clean))),
                          ")",
                          sep = "")
  
  #And now we do the same for the q:
  
  q_effect_clean <- beta_parser(tmp_end_clean$q)
  q_lower_ci_clean <- lower_ci_parser_cause(tmp_end_clean$q)
  q_upper_ci_clean <- upper_ci_parser_cause(tmp_end_clean$q)
  
  q_clean_end <- paste(as.character(round_mgu_cause(as.numeric(q_effect_clean))),
                       " (",
                       as.character(round_mgu_cause(as.numeric(q_lower_ci_clean))),
                       ", ",
                       as.character(round_mgu_cause(as.numeric(q_upper_ci_clean))),
                       ")",
                       sep = "")
  
  #And finally we do the same for the p.value
  #No need to do scientific notation with CAUSE results.
  #The lowest p-value is P=0.01.
  
  p_clean_end <- as.character(round_mgu_cause(as.numeric(tmp_end_clean$p)))
  
  tmp_end_clean_formatted <- tmp_end_clean
  
  #We need the following data for the tables...
  
  tmp_end_clean_formatted$beta <- beta_clean_end
  tmp_end_clean_formatted$q <- q_clean_end
  tmp_end_clean_formatted$p <- p_clean_end
  
  #And we need the following data for the plots
  
  tmp_end_clean_formatted$beta_effect <- beta_effect_clean
  tmp_end_clean_formatted$beta_lower_ci <- beta_lower_ci_clean
  tmp_end_clean_formatted$beta_upper_ci <- beta_upper_ci_clean
  
  #And we return:
  
  return(tmp_end_clean_formatted)
  
}

data_parsers_cause <- function(results){
  #This function is going to return a whole dataframe with the following vectors:
  #exposure
  #outcome
  #effect (either causal or shared with its confidence intervals.) This is for tables
  #beta_effect: the causal or shared effect for the plots.
  #beta_lower_ci: the lower_ci for the plots
  #beta_upper_ci: the upper_ci for the plots.
  #id (exposure-outcome)
  #pvals
  #model select
  #q = proporiton of shared variants.
  
  #First we make empty vectors to save all:
  
  outcomes <- c()
  exposures <- c()
  effects <- c()
  ids <- c()
  models <- c()
  pvals <- c()
  qs <- c()
  beta_effects <- c()
  beta_lower_cis <- c()
  beta_upper_cis <- c()
  
  #First we take the exposures!!
  
  exposure <- results$exposure_column
  exposure <- exposure[-1]
  exposure <- exposure[which(duplicated(exposure) == FALSE)]
  exposures <- rep(exposure, 5)
  
  #and we get the final results without the exposures:
  
  results <- results[,-1]
  
  #We can maybe make the column of the outcomes at this stage, so we do not get confused about it!
  
  outcome <- colnames(results)
  
  #All those that I do not need have a "."...
  #Let's use some regex to take them out...
  
  outcomes <- outcome[which(str_detect(outcome, "[.]") == FALSE)] #dots are always a bit tricky to get...
  
  #Finally,...
  
  outcomes <- rep(outcomes, 3)
  
  #With that done now we are going to clean the data.
  
  colnames_ <- results[1,]
  
  colnames(results) <- colnames_
  
  results <- results[-1,]
  
  for(i in seq(1,25, by = 5)){
    
    j=i+4
    
    tmp=results[,seq(i,j)]
    
    #NOTE WE ARE GOING TO TAKE ADVANTGE THAT ALL THE RESULTS HAVE THE SAME FORMAT.
    #WE ARE GOING TO DIVIDE THE DATA EVENLY:
    
    tmp_1 <- tmp[seq(1,2),]
    tmp_2 <- tmp[seq(3,4),]
    tmp_3 <- tmp[seq(5,6),]
    
    #We are going to clean these fellas and the joing them:
    
    tmp_1_clean <- model_parser(tmp_1)
    tmp_2_clean <- model_parser(tmp_2)
    tmp_3_clean <- model_parser(tmp_3)
    
    tmp_end_clean <- rbind(tmp_1_clean, tmp_2_clean, tmp_3_clean)
    
    ##############################################
    #From this dataframe we can already obtain...#
    ##############################################
    
    effects <- c(effects , tmp_end_clean$beta)
    beta_effects <- c(beta_effects, tmp_end_clean$beta_effect)
    beta_lower_cis <- c(beta_lower_cis, tmp_end_clean$beta_lower_ci)
    beta_upper_cis <- c(beta_upper_cis, tmp_end_clean$beta_upper_ci)
    qs <- c(qs , tmp_end_clean$q)
    pvals <- c(pvals, tmp_end_clean$p)
    models <- c(models, tmp_end_clean$model)
    #exposures <- c(exposures, c("FIadjBMI", "FI", "HOMA-IR"))
    
  }
  
  #And now we loop to get the ids
  
  for(index in seq(1, length(effects))){
    
    id <- paste(exposures[index], "-", outcomes[index], sep = "")
    
    ids <- c(ids, id)
    
  }
  
  ids <- as.data.frame(ids)
  exposures <- as.data.frame(exposures)
  outcomes <- as.data.frame(outcomes)
  effects <- as.data.frame(effects)
  qs <- as.data.frame(qs)
  pvals <- as.data.frame(pvals)
  beta_effects <- as.data.frame(beta_effects)
  beta_lower_cis <- as.data.frame(beta_lower_cis)
  beta_upper_cis <- as.data.frame(beta_upper_cis)
  models <- as.data.frame(models)
  
  final_df <- cbind(ids, exposures, outcomes, effects, beta_effects, beta_lower_cis, beta_upper_cis, qs, pvals, models)
  
  return(final_df)
  
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

mr_results_fetal <- mr_results[seq(1,24),]
mr_results_parental <- mr_results[seq(25,60),]

mr_results_fetal$Effect <- paste(unlist(sapply(as.numeric(mr_results_fetal$betas), round_mgu)), " (", unlist(sapply(as.numeric(mr_results_fetal$lower_cis), round_mgu)), ";", unlist(sapply(as.numeric(mr_results_fetal$upper_cis), round_mgu)), ")", sep = "")
mr_results_parental$Effect <- paste(unlist(sapply(as.numeric(mr_results_parental$betas), round_mgu)), " (", unlist(sapply(as.numeric(mr_results_parental$lower_cis), round_mgu)), ";", unlist(sapply(as.numeric(mr_results_parental$upper_cis), round_mgu)), ")", sep = "")

mr_results_fetal$alphabetical_id <- c("a1", "a2", "a3", "a4",
                                   "b1", "b2", "b3", "b4",
                                   "c1", "c2", "c3", "c4",
                                   "d1", "d2", "d3", "d4",
                                   "e1", "e2", "e3", "e4",
                                   "f1", "f2", "f3", "f4") 

mr_results_parental$alphabetical_id <- c("a1", "a2", "a3", "a4",
                                      "b1", "b2", "b3", "b4",
                                      "c1", "c2", "c3", "c4",
                                      "d1", "d2", "d3", "d4",
                                      "e1", "e2", "e3", "e4",
                                      "f1", "f2", "f3", "f4",
                                      "g1", "g2", "g3", "g4",
                                      "h1", "h2", "h3", "h4",
                                      "i1", "i2", "i3", "i4") 

mr_results_fetal$exposures <- factor(mr_results_fetal$exposures, levels = c("CIR", "CIRadjISI", "DI"))
mr_results_parental$exposures <- factor(mr_results_parental$exposures, levels = c("CIR", "CIRadjISI", "DI"))

mr_results_fetal$outcome_ids <- c("", "", "", "",
                             "", "", "Fetal effect on own birth weight", "",
                             "", "", "", "",
                             "", "", "", "",
                             "", "", "Fetal effect on own birth weight -adj maternal genotype", "",
                             "", "", "", "")

mr_results_parental$outcome_ids <- c("", "", "", "",
                                   "", "", "Maternal effect on offspring's birth weight", "",
                                   "", "", "", "", 
                                   "", "", "", "",
                                   "", "", "Maternal effect on offspring's birth wieght -adj fetal genotype", "",
                                   "", "", "", "",
                                   "", "", "", "",
                                   "", "", "Paternal effect on offspring's birth weight", "", 
                                   "", "", "", "")

#Now let's put the exposures, since this is gonna be easier to show:

mr_results_fetal$exposure_ids <- c("", "Fasting insulin -adj BMI", "", "",
                                  "", "Fasting insulin", "", "",
                                  "", "Insulin resistance", "", "",
                                  "", "Fasting insulin -adj BMI", "", "",
                                  "", "Fasting insulin", "", "",
                                  "", "Insulin resistance", "", "")

mr_results_parental$exposure_ids <- c("", "Fasting insulin -adj BMI", "", "",
                                     "", "Fasting insulin", "", "",
                                     "", "Insulin resistance", "", "", 
                                     "", "Fasting insulin -adj BMI", "", "",
                                     "", "Fasting insulin", "", "",
                                     "", "Insulin resistance", "", "", 
                                     "", "Fasting insulin -adj BMI", "", "",
                                     "", "Fasting insulin", "", "",
                                     "", "Insulin resistance", "", "")

#And now let's make some IDs for the methods, which is what makes this more complicated:

mr_results_fetal$methods_ids <- c("MR-Egger", "Weighted median", "IVW", "Weighted mode",
                                  "MR-Egger", "Weighted median", "IVW", "Weighted mode",
                                  "MR-Egger", "Weighted median", "IVW", "Weighted mode",
                                  "MR-Egger", "Weighted median", "IVW", "Weighted mode",
                                  "MR-Egger", "Weighted median", "IVW", "Weighted mode",
                                  "MR-Egger", "Weighted median", "IVW", "Weighted mode")

mr_results_parental$methods_ids <- c("MR-Egger", "Weighted median", "IVW", "Weighted mode",
                                     "MR-Egger", "Weighted median", "IVW", "Weighted mode",
                                     "MR-Egger", "Weighted median", "IVW", "Weighted mode",
                                     "MR-Egger", "Weighted median", "IVW", "Weighted mode",
                                     "MR-Egger", "Weighted median", "IVW", "Weighted mode",
                                     "MR-Egger", "Weighted median", "IVW", "Weighted mode",
                                     "MR-Egger", "Weighted median", "IVW", "Weighted mode",
                                     "MR-Egger", "Weighted median", "IVW", "Weighted mode",
                                     "MR-Egger", "Weighted median", "IVW", "Weighted mode")

#Let's plot the ones that are OK

mr_results_fetal <- mr_results_fetal[which(!(mr_results_fetal$lower_cis == 0 & mr_results_fetal$upper_cis == 0 & mr_results_fetal$betas == 0)),]
mr_results_parental <- mr_results_parental[which(!(mr_results_parental$lower_cis == 0 & mr_results_parental$upper_cis == 0 & mr_results_parental$betas == 0)),]

mr_plot_fetal <- ggplot(mr_results_fetal, aes(y = fct_rev(alphabetical_id), x = as.numeric(betas), colour = exposures)) +
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
        legend.position = "none")


mr_plot_parental <- ggplot(mr_results_parental, aes(y = fct_rev(alphabetical_id), x = as.numeric(betas), colour = exposures)) +
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
        legend.position = "none")



table_base_fetal <- ggplot(mr_results_fetal, aes(y=alphabetical_id)) +
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

table_base_parental <- ggplot(mr_results_parental, aes(y=alphabetical_id)) +
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

#Let's make some tables:

tab_outcome_fetal <- table_base_fetal + 
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = outcome_ids), size = 3, lineheight = .2,
            inherit.aes = TRUE, vjust = 1)

tab_outcome_parental <- table_base_parental + 
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = outcome_ids), size = 3, lineheight = .2,
            inherit.aes = TRUE, vjust = 1)

#Now the table with the exposures:

tab_exposure_fetal <- table_base_fetal + 
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = exposure_ids), size = 3, lineheight = .2,
            inherit.aes = TRUE, vjust = 1)

tab_exposure_parental <- table_base_parental + 
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = exposure_ids), size = 3, lineheight = .2,
            inherit.aes = TRUE, vjust = 1)

#And now with the methods:

tab_methods_fetal <- table_base_fetal + 
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = methods_ids), size = 3, lineheight = .2,
            inherit.aes = TRUE, vjust = 1)

tab_methods_parental <- table_base_parental + 
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = methods_ids), size = 3, lineheight = .2,
            inherit.aes = TRUE, vjust = 1)


## 95% CI table
tab1_fetal <- table_base_fetal +
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = Effect), size = 2.7,lineheight = 1,
            inherit.aes = TRUE,vjust = 1)

## 95% CI table
tab1_parental <- table_base_parental +
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = Effect), size = 2.7,lineheight = 1,
            inherit.aes = TRUE,vjust = 1)

#pval tables
tab2_fetal <- table_base_fetal +
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = unlist(sapply(as.numeric(pvals), round_mgu_pval))), size = 2.7,lineheight = 1,
            inherit.aes = TRUE, vjust = 1)

#pval tables
tab2_parental <- table_base_parental +
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = unlist(sapply(as.numeric(pvals), round_mgu_pval))), size = 2.7,lineheight = 1,
            inherit.aes = TRUE, vjust = 1)

lay_fetal <-  matrix(c(1,1,1,2,2,2,3,3,3,3,3,3,3,4,4,4,5,5,5), nrow = 1)
lay_parental <-  matrix(c(1,1,1,2,2,2,3,3,3,3,3,3,3,4,4,4,5,5,5), nrow = 1)

plot_fetal <- grid.arrange(tab_outcome_fetal, tab_methods_fetal, mr_plot_fetal,  tab1_fetal, tab2_fetal, layout_matrix = lay_fetal)
plot_parental <- grid.arrange(tab_outcome_parental, tab_methods_parental, mr_plot_parental, tab1_parental, tab2_parental, layout_matrix = lay_parental)

#Though we are gonna add CAUSE now, it is better if we print it at least once:

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/figures/2SMR/IS_2SMR_plot.tiff", width = 6000, height=8000, res = 300)
grid.arrange(plot_fetal, plot_parental)
dev.off()

##############################################################################################################
##############################################################################################################
##############################################################################################################

cause_results <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/supplementary_tables_per_method/CAUSE/raw/CAUSE_IS.csv")

cause_results_og <- as.data.frame(cause_results)

cause_results <- data_parsers_cause(cause_results_og)

#We are going to try to add the CAUSE data to the mr_dataframes we had before
#This might allow us to make plots that can be combined together.
#Once that is done, we just need to report them in the results methods.

#First let's add the beta effects:

mr_results_fetal$cause_betas_effects <- c(cause_results$beta_effects[1],
                                  3, 3, cause_results$beta_effects[2], 3,
                                  3, 3, cause_results$beta_effects[3], 3,
                                  3, 3, cause_results$beta_effects[4], 3,
                                  cause_results$beta_effects[5], cause_results$beta_effects[6])

mr_results_parental$cause_betas_effects <- c(3, 3, cause_results$beta_effects[7], 3,
                                     3, 3, cause_results$beta_effects[8], 3,
                                     3, 3, cause_results$beta_effects[9], 3, 
                                     3, 3, cause_results$beta_effects[10], 3,
                                     3, 3, cause_results$beta_effects[11], 3,
                                     3, 3, cause_results$beta_effects[12], 3,
                                     3, 3, cause_results$beta_effects[13], 3,
                                     3, 3, cause_results$beta_effects[14], 3, 
                                     3, 3, cause_results$beta_effects[15], 3)

#Now we add the lower and upper cis:

mr_results_fetal$cause_betas_lower_cis <- c(cause_results$beta_lower_cis[1],
                                          3, 3, cause_results$beta_lower_cis[2], 3,
                                          3, 3, cause_results$beta_lower_cis[3], 3,
                                          3, 3, cause_results$beta_lower_cis[4], 3,
                                          cause_results$beta_lower_cis[5], cause_results$beta_lower_cis[6])

mr_results_fetal$cause_betas_upper_cis <- c(cause_results$beta_upper_cis[1],
                                            3, 3, cause_results$beta_upper_cis[2], 3,
                                            3, 3, cause_results$beta_upper_cis[3], 3,
                                            3, 3, cause_results$beta_upper_cis[4], 3,
                                            cause_results$beta_upper_cis[5], cause_results$beta_upper_cis[6])

mr_results_parental$cause_betas_lower_cis <- c(3, 3, cause_results$beta_lower_cis[7], 3,
                                            3, 3, cause_results$beta_lower_cis[8], 3,
                                            3, 3, cause_results$beta_lower_cis[9], 3, 
                                            3, 3, cause_results$beta_lower_cis[10], 3,
                                            3, 3, cause_results$beta_lower_cis[11], 3,
                                            3, 3, cause_results$beta_lower_cis[12], 3,
                                            3, 3, cause_results$beta_lower_cis[13], 3,
                                            3, 3, cause_results$beta_lower_cis[14], 3, 
                                            3, 3, cause_results$beta_lower_cis[15], 3)

mr_results_parental$cause_betas_upper_cis <- c(3, 3, cause_results$beta_upper_cis[7], 3,
                                              3, 3, cause_results$beta_upper_cis[8], 3,
                                              3, 3, cause_results$beta_upper_cis[9], 3, 
                                              3, 3, cause_results$beta_upper_cis[10], 3,
                                              3, 3, cause_results$beta_upper_cis[11], 3,
                                              3, 3, cause_results$beta_upper_cis[12], 3,
                                              3, 3, cause_results$beta_upper_cis[13], 3,
                                              3, 3, cause_results$beta_upper_cis[14], 3, 
                                              3, 3, cause_results$beta_upper_cis[15], 3)

#Now let's get the effects that will go to the tables:

mr_results_fetal$cause_effects <- c(cause_results$effects[1],
                                          "", "", cause_results$effects[2], "",
                                          "", "", cause_results$effects[3], "",
                                          "", "", cause_results$effects[4], "",
                                          cause_results$effects[5], cause_results$effects[6])

mr_results_parental$cause_effects <- c("", "", cause_results$effects[7], "",
                                            "", "", cause_results$effects[8], "",
                                            "", "", cause_results$effects[9], "", 
                                            "", "", cause_results$effects[10], "",
                                            "", "", cause_results$effects[11], "",
                                            "", "", cause_results$effects[12], "",
                                            "", "", cause_results$effects[13], "",
                                            "", "", cause_results$effects[14], "", 
                                            "", "", cause_results$effects[15], "")

#We are almost there, now the pvals:

mr_results_fetal$cause_pvals <- c(cause_results$pvals[1],
                                    "", "", cause_results$pvals[2], "",
                                    "", "", cause_results$pvals[3], "",
                                    "", "", cause_results$pvals[4], "",
                                    cause_results$pvals[5], cause_results$pvals[6])

mr_results_parental$cause_pvals <- c("", "", cause_results$pvals[7], "",
                                       "", "", cause_results$pvals[8], "",
                                       "", "", cause_results$pvals[9], "", 
                                       "", "", cause_results$pvals[10], "",
                                       "", "", cause_results$pvals[11], "",
                                       "", "", cause_results$pvals[12], "",
                                       "", "", cause_results$pvals[13], "",
                                       "", "", cause_results$pvals[14], "", 
                                       "", "", cause_results$pvals[15], "")

#Finally, the most important of all: the method selected:

mr_results_fetal$cause_models <- c(cause_results$models[1],
                                  "", "", cause_results$models[2], "",
                                  "", "", cause_results$models[3], "",
                                  "", "", cause_results$models[4], "",
                                  cause_results$models[5], cause_results$models[6])

mr_results_parental$cause_models <- c("", "", cause_results$models[7], "",
                                     "", "", cause_results$models[8], "",
                                     "", "", cause_results$models[9], "", 
                                     "", "", cause_results$models[10], "",
                                     "", "", cause_results$models[11], "",
                                     "", "", cause_results$models[12], "",
                                     "", "", cause_results$models[13], "",
                                     "", "", cause_results$models[14], "", 
                                     "", "", cause_results$models[15], "")


#########################################
#NOW WE CAN GENERATE THE TWO CAUSE PLOTS#
#########################################

cause_plot_fetal <- ggplot(mr_results_fetal, aes(y = fct_rev(alphabetical_id), x = as.numeric(cause_betas_effects), colour = exposures)) +
  geom_point(size = 2.7, position = "identity") +  
  geom_errorbarh(aes(xmin = as.numeric(cause_betas_lower_cis), xmax = as.numeric(cause_betas_upper_cis)), height = 0.25) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  xlim(-2, 2)+
  xlab("Median causal/shared effect (95% CI)") + 
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


cause_plot_parental <- ggplot(mr_results_parental, aes(y = fct_rev(alphabetical_id), x = as.numeric(cause_betas_effects), colour = exposures)) +
  geom_point(size = 2.7, position = "identity") +  
  geom_errorbarh(aes(xmin = as.numeric(cause_betas_lower_cis), xmax = as.numeric(cause_betas_upper_cis)), height = 0.25) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  xlim(-2.99, 2.99)+
  xlab("Median causal/shared effect (95% CI)") + 
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

#####################
#Let's do the tables#
#####################

table_cause_fetal <- ggplot(mr_results_fetal, aes(y=alphabetical_id)) +
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

table_cause_parental <- ggplot(mr_results_parental, aes(y=alphabetical_id)) +
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


#And now with the methods:

tab_cause_model_fetal <- table_cause_fetal + 
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = cause_models), size = 3, lineheight = .2,
            inherit.aes = TRUE, vjust = 1)

tab_cause_model_parental <- table_cause_parental + 
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = cause_models), size = 3, lineheight = .2,
            inherit.aes = TRUE, vjust = 1)


## 95% CI table
tab_cause_effects_fetal <- table_cause_fetal +
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = cause_effects), size = 2.7,lineheight = 1,
            inherit.aes = TRUE,vjust = 1)

## 95% CI table
tab_cause_effects_parental <- table_cause_parental +
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = cause_effects), size = 2.7,lineheight = 1,
            inherit.aes = TRUE,vjust = 1)

#pval tables
tab_cause_pval_fetal <- table_cause_fetal +
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = unlist(sapply(as.numeric(cause_pvals), round_mgu_pval))), size = 2.7,lineheight = 1,
            inherit.aes = TRUE, vjust = 1)

#pval tables
tab_cause_pval_parental <- table_cause_parental +
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = unlist(sapply(as.numeric(cause_pvals), round_mgu_pval))), size = 2.7,lineheight = 1,
            inherit.aes = TRUE, vjust = 1)

lay_fetal <-  matrix(c(1,1,1,1,1,1,1,1,2,2,2,3,3,3,3,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,8,8,8,8,8,9,9,9), nrow = 1)
lay_parental <-  matrix(c(1,1,1,1,1,1,1,1,1,2,2,2,3,3,3,3,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,8,8,8,8,8,9,9,9), nrow = 1)

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/figures/2SMR/IS_fetal_2SMR_CAUSE_plot.tiff", width = 6000, height=1000, res = 300)
plot_fetal <- grid.arrange(tab_outcome_fetal, tab_methods_fetal, mr_plot_fetal,  tab1_fetal, tab2_fetal, tab_cause_model_fetal, tab_cause_pval_fetal, cause_plot_fetal, tab_cause_effects_fetal, layout_matrix = lay_fetal)
dev.off()

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/figures/2SMR/IS_parental_2SMR_CAUSE_plot.tiff", width = 6000, height=2000, res = 300)
plot_parental <- grid.arrange(tab_outcome_parental, tab_methods_parental, mr_plot_parental, tab1_parental, tab2_parental, tab_cause_model_parental, tab_cause_pval_parental, cause_plot_parental, tab_cause_effects_parental, layout_matrix = lay_parental)
dev.off()

