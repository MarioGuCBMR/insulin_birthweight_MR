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

###################
#Loading functions#
###################

lower_ci <- function(beta_, se_){
  
  lower_ci <- beta_ - qnorm(0.975)*se_
  
  return(lower_ci)
  
}

upper_ci <- function(beta_, se_){
  
  upper_ci <- beta_ + qnorm(0.975)*se_
  
  return(upper_ci)
  
}

##############
#Loading data#
##############

CIR_fbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/CIR/FBW/CIR_FBW_weighted")
CIR_mbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/CIR/MBW/CIR_MBW_weighted")
CIR_fbwadjmbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/CIR/FBWadjMBW/CIR_FBWadjMBW_weighted")
CIR_mbwadjfbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/CIR/MBWadjFBW/CIR_mbwadjfbw_weighted")
CIR_father <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/CIR/Father/CIR_father_weighted")

#CIRadjISI:

CIRadjISI_fbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/CIRadjISI/FBW/CIRadjISI_FBW_weighted")
CIRadjISI_mbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/CIRadjISI/MBW/CIRadjISI_MBW_weighted")
CIRadjISI_fbwadjmbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/CIRadjISI/FBWadjMBW/CIRadjISI_FBWadjMBW_weighted")
CIRadjISI_mbwadjfbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/CIRadjISI/MBWadjFBW/CIRadjISI_mbwadjfbw_weighted")
CIRadjISI_father <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/CIRadjISI/Father/CIRadjISI_father_weighted")

#DI:

DI_fbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/DI/FBW/DI_FBW_weighted")
DI_mbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/DI/MBW/DI_MBW_weighted")
DI_fbwadjmbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/DI/FBWadjMBW/DI_FBWadjMBW_weighted")
DI_mbwadjfbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/DI/MBWadjFBW/DI_mbwadjfbw_weighted")
DI_father <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/DI/Father/DI_father_weighted")

##################################################
#Let's make a dataframe with the betas and the SE#
##################################################

betas_ <- c(CIR_fbw$ahat, CIRadjISI_fbw$ahat, DI_fbw$ahat, CIR_fbwadjmbw$ahat, CIRadjISI_fbwadjmbw$ahat, DI_fbwadjmbw$ahat,
           CIR_mbw$ahat, CIRadjISI_mbw$ahat, DI_mbw$ahat, CIR_mbwadjfbw$ahat, CIRadjISI_mbwadjfbw$ahat, DI_mbwadjfbw$ahat,
           CIR_father$ahat, CIRadjISI_father$ahat, DI_father$ahat)

se_ <- c(CIR_fbw$aSE, CIRadjISI_fbw$aSE, DI_fbw$aSE, CIR_fbwadjmbw$aSE, CIRadjISI_fbwadjmbw$aSE, DI_fbwadjmbw$aSE,
           CIR_mbw$aSE, CIRadjISI_mbw$aSE, DI_mbw$aSE, CIR_mbwadjfbw$aSE, CIRadjISI_mbwadjfbw$aSE, DI_mbwadjfbw$aSE,
           CIR_father$aSE, CIRadjISI_father$aSE, DI_father$aSE)

pval_ <-  c(CIR_fbw$pval, CIRadjISI_fbw$pval, DI_fbw$pval, CIR_fbwadjmbw$pval, CIRadjISI_fbwadjmbw$pval, DI_fbwadjmbw$pval,
            CIR_mbw$pval, CIRadjISI_mbw$pval, DI_mbw$pval, CIR_mbwadjfbw$pval, CIRadjISI_mbwadjfbw$pval, DI_mbwadjfbw$pval,
            CIR_father$pval, CIRadjISI_father$pval, DI_father$pval)

exposure_ <-  c("CIR", "CIRadjISI", "DI", "CIR", "CIRadjISI", "DI",
            "CIR", "CIRadjISI", "DI", "CIR", "CIRadjISI", "DI",
            "CIR", "CIRadjISI", "DI")

outcome_ <-  c("Fetal effects on own birthweight", "Fetal effects on own birthweight", "Fetal effects on own birthweight",
               "Fetal effects on own birthweight adjusted for maternal effect", "Fetal effects on own birthweight adjusted for maternal effect", "Fetal effects on own birthweight adjusted for maternal effect",
               "Maternal effects on offspring's birthweight", "Maternal effects on offspring's birthweight", "Maternal effects on offspring's birthweight",
               "Maternal effects on offspring's birthweight adjusted for fetal effect", "Maternal effects on offspring's birthweight adjusted for fetal effect", "Maternal effects on offspring's birthweight adjusted for fetal effect",
               "Paternal effects on offspring's birthweight", "Paternal effects on offspring's birthweight", "Paternal effects on offspring's birthweight")



df <- cbind(as.data.frame(betas_), as.data.frame(se_), as.data.frame(pval_), as.data.frame(exposure_), as.data.frame(outcome_))

colnames(df) <- c("Beta", "SE", "P", "Exposure", "Outcome")

df$lower_ci <- "-"
df$upper_ci <- "-"

for(i in seq(1, length(df$Beta))){
  
  df$lower_ci[i] <- lower_ci(df$Beta[i], df$SE[i])
  df$upper_ci[i] <- upper_ci(df$Beta[i], df$SE[i])
  
}

df$effect <- paste(round(as.numeric(df$Beta), digits = 2), " (", round(as.numeric(df$lower_ci), digits = 2), ";", round(as.numeric(df$upper_ci), digits = 2), ")", sep = "")
df$P <- formatC(df$P, format = "e", digits = 2)

############
#Make plots#
############

Outcome_order <- outcome_

dotCOLS = c("red3","red1", "deepskyblue3", "deepskyblue1", "green3")
barCOLS = c("red4","red2", "deepskyblue4", "deepskyblue2", "green4")

p <- ggplot(df, aes(x=fct_rev(Outcome), y=as.numeric(Beta), ymin=as.numeric(lower_ci), ymax=as.numeric(upper_ci),col=Exposure,fill=Exposure)) + 
  #specify position here
  geom_errorbar(size=0.5,position=position_dodge(width = -.9)) +
  geom_hline(yintercept=0) +
  #specify position here too
  geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = -.9)) +
  scale_fill_manual(values=barCOLS)+
  scale_color_manual(values=dotCOLS)+
  scale_x_discrete(name="Outcome") +
  scale_y_continuous(name = "GRS effect size with 95% confidence intervals", breaks = c(-0.20, -0.10, 0, 0.10, 0.20), limits = c(-0.25, 0.25)) +
  coord_flip() +
  theme_minimal() +
  theme(legend.position="bottom")

df_table <- df %>%
  select(Outcome, effect, P)

df_table$Outcome <- seq(1,8, by = 0.5)

colour <- "white"

df_table$colour <- colour

data_table <- ggplot(data=df_table, aes(y= as.numeric(fct_rev(as.character(Outcome)))))+
  geom_text(aes(x= 1, label=effect), hjust = 1, size = 3)+
  geom_text(aes(x= 1.2, label=P), hjust = 0, size = 3)+
  scale_x_continuous(name = "", breaks = c(0.5, 1.5), limits = c(0.5, 1.5))+
  scale_y_continuous(name = "", breaks = seq(1,15, by = 0.5), limits = c(1, 15))+
  theme(axis.text.y = element_text(lineheight = 0.5, size = 6)) +
  theme_void()
  

pdp::grid.arrange(p, data_table, ncol = 2, widths = c(2.5,1))
