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
#LoaIR_53ng functions#
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
#LoaIR_53ng data#
##############

FIadjBMI_fbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/FIadjBMI/FBW/FIadjBMI_FBW_weighted")
FIadjBMI_mbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/FIadjBMI/MBW/FIadjBMI_MBW_weighted")
FIadjBMI_fbwadjmbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/FIadjBMI/FBWadjMBW/FIadjBMI_FBWadjMBW_weighted")
FIadjBMI_mbwadjfbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/FIadjBMI/MBWadjFBW/FIadjBMI_mbwadjfbw_weighted")
FIadjBMI_father <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/FIadjBMI/Father/FIadjBMI_father_weighted")

#FI:

#FI_fbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/FI/FBW/FI_FBW_weighted")
#FI_mbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/FI/MBW/FI_MBW_weighted")
#FI_fbwadjmbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/FI/FBWadjMBW/FI_FBWadjMBW_weighted")
#FI_mbwadjfbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/FI/MBWadjFBW/FI_mbwadjfbw_weighted")
#FI_father <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/FI/Father/FI_father_weighted")

#IR_53:

IR_53_fbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/53_IR_Lotta_SNPs/FBW/IR_53_FBW_weighted")
IR_53_mbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/53_IR_Lotta_SNPs/MBW/IR_53_MBW_weighted")
IR_53_fbwadjmbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/53_IR_Lotta_SNPs/FBWadjMBW/IR_53_FBWadjMBW_weighted")
IR_53_mbwadjfbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/53_IR_Lotta_SNPs/MBWadjFBW/IR_53_mbwadjfbw_weighted")
IR_53_father <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/GRS/53_IR_Lotta_SNPs/Father/IR_53_father_weighted")

##################################################
#Let's make a dataframe with the betas and the SE#
##################################################

betas_ <- c(FIadjBMI_fbw$ahat, IR_53_fbw$ahat, FIadjBMI_fbwadjmbw$ahat, IR_53_fbwadjmbw$ahat,
           FIadjBMI_mbw$ahat, IR_53_mbw$ahat, FIadjBMI_mbwadjfbw$ahat, IR_53_mbwadjfbw$ahat,
           FIadjBMI_father$ahat, IR_53_father$ahat)

se_ <- c(FIadjBMI_fbw$aSE, IR_53_fbw$aSE, FIadjBMI_fbwadjmbw$aSE, IR_53_fbwadjmbw$aSE,
           FIadjBMI_mbw$aSE, IR_53_mbw$aSE, FIadjBMI_mbwadjfbw$aSE, IR_53_mbwadjfbw$aSE,
           FIadjBMI_father$aSE, IR_53_father$aSE)

pval_ <-  c(FIadjBMI_fbw$pval, IR_53_fbw$pval, FIadjBMI_fbwadjmbw$pval, IR_53_fbwadjmbw$pval,
            FIadjBMI_mbw$pval, IR_53_mbw$pval, FIadjBMI_mbwadjfbw$pval, IR_53_mbwadjfbw$pval,
            FIadjBMI_father$pval, IR_53_father$pval)

exposure_ <-  c("FIadjBMI", "IR_53", "FIadjBMI","IR_53",
            "FIadjBMI", "IR_53", "FIadjBMI", "IR_53",
            "FIadjBMI", "IR_53")

outcome_ <-  c("Fetal effects on own birthweight", "Fetal effects on own birthweight",
               "Fetal effects on own birthweight adjusted for maternal effect", "Fetal effects on own birthweight adjusted for maternal effect",
               "Maternal effects on offspring's birthweight", "Maternal effects on offspring's birthweight",
               "Maternal effects on offspring's birthweight adjusted for fetal effect", "Maternal effects on offspring's birthweight adjusted for fetal effect",
               "Paternal effects on offspring's birthweight", "Paternal effects on offspring's birthweight")



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
  scale_y_continuous(name = "GRS effect size with 95% confidence intervals", breaks = c(-0.50, -0.40, -0.30, -0.20, -0.10, 0, 0.10, 0.20, 0.30, 0.40, 0.50), limits = c(-0.50, 0.50)) +
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
