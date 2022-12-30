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

beta_parser <- function(effect_size){
  
  tmp <- as.character(strsplit(effect_size, " ")[[1]][1])
  
  return(tmp)
  
}

lower_ci_parser <- function(effect_size){
  
  tmp <- as.character(strsplit(effect_size, " ")[[1]][2])
  tmp <- as.character(strsplit(tmp, "[()]")[[1]][2])
  tmp <- as.character(strsplit(tmp, ",")[[1]][1])
  
  
  return(tmp)
  
}

upper_ci_parser <- function(effect_size){
  
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
  beta_lower_ci_clean <- lower_ci_parser(tmp_end_clean$beta)
  beta_upper_ci_clean <- upper_ci_parser(tmp_end_clean$beta)
  
  beta_clean_end <- paste(as.character(round_mgu(as.numeric(beta_effect_clean))),
                          " (",
                          as.character(round_mgu(as.numeric(beta_lower_ci_clean))),
                          ", ",
                          as.character(round_mgu(as.numeric(beta_upper_ci_clean))),
                          ")",
                          sep = "")
    
  #And now we do the same for the q:
  
  q_effect_clean <- beta_parser(tmp_end_clean$q)
  q_lower_ci_clean <- lower_ci_parser(tmp_end_clean$q)
  q_upper_ci_clean <- upper_ci_parser(tmp_end_clean$q)
  
  q_clean_end <- paste(as.character(round_mgu(as.numeric(q_effect_clean))),
                          " (",
                          as.character(round_mgu(as.numeric(q_lower_ci_clean))),
                          ", ",
                          as.character(round_mgu(as.numeric(q_upper_ci_clean))),
                          ")",
                          sep = "")
  
  #And finally we do the same for the p.value
  #No need to do scientific notation with CAUSE results.
  #The lowest p-value is P=0.01.
  
  p_clean_end <- as.character(round_mgu(as.numeric(tmp_end_clean$p)))
  
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

data_parsers <- function(results){
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

cause_results <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/supplementary_tables_per_method/CAUSE/raw/CAUSE_IR.csv")

cause_results_og <- as.data.frame(cause_results)

cause_results <- data_parsers(cause_results_og)

#################################
#Make plots for weighted section#
#################################

cause_results$alphabetical_id <- c("a1", "a2", "a3", "a4", "a5",
                                         "b1", "b2", "b3", "b4", "b5",
                                         "c1", "c2", "c3", "c4", "c5")

cause_results$exposures <- factor(cause_results$exposures, levels = c("FIadjBMI", "FI", "HOMA-IR"))

plot1 <- ggplot(cause_results, aes(y = fct_rev(alphabetical_id), x = as.numeric(beta_effects), shape = exposures)) +
  geom_point(size = 2.7, position = "identity") +  
  geom_errorbarh(aes(xmin = as.numeric(beta_lower_cis), xmax = as.numeric(beta_upper_cis)), height = 0.25) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  xlab("Effect size (95% CI)") + 
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


cause_results$final_ids <- c("", "FBW", "",
                                   "", "FBWadjMBW", "",
                                   "", "MBW", "",
                                   "", "MBWadjFBW", "",
                                   "", "Father", "")

table_base <- ggplot(cause_results, aes(y=alphabetical_id)) +
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

#Outcome table

tab0 <- table_base + 
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = final_ids), size = 3, lineheight = .2,
            inherit.aes = TRUE, vjust = 1)

#Model tables

tab1 <- table_base + 
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = models), size = 3, lineheight = .2,
            inherit.aes = TRUE, vjust = 1)

#Pvals

tab2 <- table_base + 
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = pvals), size =2.7,lineheight = 1,
            inherit.aes = TRUE, vjust = 1)


## q values (maybe we will remove them...)
tab3 <- table_base +
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = qs), size = 2.7,lineheight = 1,
            inherit.aes = TRUE,vjust = 1)

#the real effects sizes
tab4 <- table_base +
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = effects), size = 2.7,lineheight = 1,
            inherit.aes = TRUE, vjust = 1)


lay <-  matrix(c(1,1,1,2,2,3,3,3,4,4,4,5,5,5,6,6,6,6,6,6,6,6,6,6), nrow = 1)

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/figures/CAUSE/IR_CAUSE_plot.tiff", width = 3000, height=1200, res = 300)
grid.arrange(tab0, tab1, tab2, tab3, tab4, plot1, layout_matrix = lay)
dev.off()
