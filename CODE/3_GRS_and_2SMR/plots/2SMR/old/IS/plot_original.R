##############
#INTRODUCTION#
##############

#This code will take all the 2SMR results and make forest plots automatically.
#This is going to take the original 2SMR only.
#We will do the unoriginal later.

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

CIR_fbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/CIR/FBW/CIR_FBW_original")
CIR_mbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/CIR/MBW/CIR_MBW_original")
CIR_fbwadjmbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/CIR/FBWadjMBW/CIR_FBWadjMBW_original")
CIR_mbwadjfbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/CIR/MBWadjFBW/CIR_mbwadjfbw_original")
CIR_father <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/CIR/Father/CIR_father_original")

#CIRadjISI:

CIRadjISI_fbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/CIRadjISI/FBW/CIRadjISI_FBW_original")
CIRadjISI_mbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/CIRadjISI/MBW/CIRadjISI_MBW_original")
CIRadjISI_fbwadjmbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/CIRadjISI/FBWadjMBW/CIRadjISI_FBWadjMBW_original")
CIRadjISI_mbwadjfbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/CIRadjISI/MBWadjFBW/CIRadjISI_mbwadjfbw_original")
CIRadjISI_father <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/CIRadjISI/Father/CIRadjISI_father_original")

#DI:

DI_fbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/DI/FBW/DI_FBW_original")
DI_mbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/DI/MBW/DI_MBW_original")
DI_fbwadjmbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/DI/FBWadjMBW/DI_FBWadjMBW_original")
DI_mbwadjfbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/DI/MBWadjFBW/DI_mbwadjfbw_original")
DI_father <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/DI/Father/DI_father_original")

#################################################
#Let's clean simple mode because nobody likes it#
#################################################

CIR_fbw <- CIR_fbw[which(CIR_fbw$method != "Simple mode"),]
CIR_fbwadjmbw <- CIR_fbwadjmbw[which(CIR_fbwadjmbw$method != "Simple mode"),]
CIR_mbw <- CIR_mbw[which(CIR_mbw$method != "Simple mode"),]
CIR_mbwadjfbw <- CIR_mbwadjfbw[which(CIR_mbwadjfbw$method != "Simple mode"),]
CIR_father <- CIR_father[which(CIR_father$method != "Simple mode"),]

CIRadjISI_fbw <- CIRadjISI_fbw[which(CIRadjISI_fbw$method != "Simple mode"),]
CIRadjISI_fbwadjmbw <- CIRadjISI_fbwadjmbw[which(CIRadjISI_fbwadjmbw$method != "Simple mode"),]
CIRadjISI_mbw <- CIRadjISI_mbw[which(CIRadjISI_mbw$method != "Simple mode"),]
CIRadjISI_mbwadjfbw <- CIRadjISI_mbwadjfbw[which(CIRadjISI_mbwadjfbw$method != "Simple mode"),]
CIRadjISI_father <- CIRadjISI_father[which(CIRadjISI_father$method != "Simple mode"),]

DI_fbw <- DI_fbw[which(DI_fbw$method != "Simple mode"),]
DI_fbwadjmbw <- DI_fbwadjmbw[which(DI_fbwadjmbw$method != "Simple mode"),]
DI_mbw <- DI_mbw[which(DI_mbw$method != "Simple mode"),]
DI_mbwadjfbw <- DI_mbwadjfbw[which(DI_mbwadjfbw$method != "Simple mode"),]
DI_father <- DI_father[which(DI_father$method != "Simple mode"),]

##################################################
#Let's make a dataframe with the betas and the SE#
##################################################

betas_ <- c(CIR_fbw$b, CIRadjISI_fbw$b, DI_fbw$b, CIR_fbwadjmbw$b, CIRadjISI_fbwadjmbw$b, DI_fbwadjmbw$b,
           CIR_mbw$b, CIRadjISI_mbw$b, DI_mbw$b, CIR_mbwadjfbw$b, CIRadjISI_mbwadjfbw$b, DI_mbwadjfbw$b,
           CIR_father$b, CIRadjISI_father$b, DI_father$b)

se_ <- c(CIR_fbw$se, CIRadjISI_fbw$se, DI_fbw$se, CIR_fbwadjmbw$se, CIRadjISI_fbwadjmbw$se, DI_fbwadjmbw$se,
           CIR_mbw$se, CIRadjISI_mbw$se, DI_mbw$se, CIR_mbwadjfbw$se, CIRadjISI_mbwadjfbw$se, DI_mbwadjfbw$se,
           CIR_father$se, CIRadjISI_father$se, DI_father$se)

pval_ <-  c(CIR_fbw$pval, CIRadjISI_fbw$pval, DI_fbw$pval, CIR_fbwadjmbw$pval, CIRadjISI_fbwadjmbw$pval, DI_fbwadjmbw$pval,
            CIR_mbw$pval, CIRadjISI_mbw$pval, DI_mbw$pval, CIR_mbwadjfbw$pval, CIRadjISI_mbwadjfbw$pval, DI_mbwadjfbw$pval,
            CIR_father$pval, CIRadjISI_father$pval, DI_father$pval)

exposure_ <-  c("CIR", "CIR", "CIR", "CIR", "CIRadjISI", "CIRadjISI", "CIRadjISI", "CIRadjISI", "DI", "DI", "DI", "DI",
                "CIR", "CIR", "CIR", "CIR", "CIRadjISI", "CIRadjISI", "CIRadjISI", "CIRadjISI", "DI", "DI", "DI", "DI",
                "CIR", "CIR", "CIR", "CIR", "CIRadjISI", "CIRadjISI", "CIRadjISI", "CIRadjISI", "DI", "DI", "DI", "DI",
                "CIR", "CIR", "CIR", "CIR", "CIRadjISI", "CIRadjISI", "CIRadjISI", "CIRadjISI", "DI", "DI", "DI", "DI",
                "CIR", "CIR", "CIR", "CIR", "CIRadjISI", "CIRadjISI", "CIRadjISI", "CIRadjISI", "DI", "DI", "DI", "DI")

method_ <-  c("Egger-CIR", "Weigthed_Median_CIR", "IVW_CIR", "Weighted_Mode_CIR", "Egger-CIRadjISI", "Weighted_Median_CIRadjISI", "Egger_CIRadjISI", "Weighted_Mode_CIRadjISI", "Egger_DI", "Weighted_Median_DI", "IVW_DI", "Weighted_Mode_DI",
              "Egger-CIR", "Weigthed_Median_CIR", "IVW_CIR", "Weighted_Mode_CIR", "Egger-CIRadjISI", "Weighted_Median_CIRadjISI", "Egger_CIRadjISI", "Weighted_Mode_CIRadjISI", "Egger_DI", "Weighted_Median_DI", "IVW_DI", "Weighted_Mode_DI",
              "Egger-CIR", "Weigthed_Median_CIR", "IVW_CIR", "Weighted_Mode_CIR", "Egger-CIRadjISI", "Weighted_Median_CIRadjISI", "Egger_CIRadjISI", "Weighted_Mode_CIRadjISI", "Egger_DI", "Weighted_Median_DI", "IVW_DI", "Weighted_Mode_DI",
              "Egger-CIR", "Weigthed_Median_CIR", "IVW_CIR", "Weighted_Mode_CIR", "Egger-CIRadjISI", "Weighted_Median_CIRadjISI", "Egger_CIRadjISI", "Weighted_Mode_CIRadjISI", "Egger_DI", "Weighted_Median_DI", "IVW_DI", "Weighted_Mode_DI",
              "Egger-CIR", "Weigthed_Median_CIR", "IVW_CIR", "Weighted_Mode_CIR", "Egger-CIRadjISI", "Weighted_Median_CIRadjISI", "Egger_CIRadjISI", "Weighted_Mode_CIRadjISI", "Egger_DI", "Weighted_Median_DI", "IVW_DI", "Weighted_Mode_DI")

MR <-  c(CIR_fbw$method, CIR_fbw$method, CIR_fbw$method, CIR_fbw$method, CIR_fbw$method,
         CIR_fbw$method, CIR_fbw$method, CIR_fbw$method, CIR_fbw$method, CIR_fbw$method,
         CIR_fbw$method, CIR_fbw$method, CIR_fbw$method, CIR_fbw$method, CIR_fbw$method)


outcome_ <-  c("Fetal effects on own birthweight", "Fetal effects on own birthweight", "Fetal effects on own birthweight", "Fetal effects on own birthweight",
               "Fetal effects on own birthweight", "Fetal effects on own birthweight", "Fetal effects on own birthweight", "Fetal effects on own birthweight",
               "Fetal effects on own birthweight", "Fetal effects on own birthweight", "Fetal effects on own birthweight", "Fetal effects on own birthweight",
               "Fetal effects on own birthweight adjusted for maternal effect", "Fetal effects on own birthweight adjusted for maternal effect", "Fetal effects on own birthweight adjusted for maternal effect", "Fetal effects on own birthweight adjusted for maternal effect",
               "Fetal effects on own birthweight adjusted for maternal effect", "Fetal effects on own birthweight adjusted for maternal effect", "Fetal effects on own birthweight adjusted for maternal effect", "Fetal effects on own birthweight adjusted for maternal effect",
               "Fetal effects on own birthweight adjusted for maternal effect", "Fetal effects on own birthweight adjusted for maternal effect", "Fetal effects on own birthweight adjusted for maternal effect", "Fetal effects on own birthweight adjusted for maternal effect",
               "Maternal effects on offspring's birthweight", "Maternal effects on offspring's birthweight", "Maternal effects on offspring's birthweight", "Maternal effects on offspring's birthweight",
               "Maternal effects on offspring's birthweight", "Maternal effects on offspring's birthweight", "Maternal effects on offspring's birthweight", "Maternal effects on offspring's birthweight",
               "Maternal effects on offspring's birthweight", "Maternal effects on offspring's birthweight", "Maternal effects on offspring's birthweight", "Maternal effects on offspring's birthweight",
               "Maternal effects on offspring's birthweight adjusted for fetal effect", "Maternal effects on offspring's birthweight adjusted for fetal effect", "Maternal effects on offspring's birthweight adjusted for fetal effect", "Maternal effects on offspring's birthweight adjusted for fetal effect",
               "Maternal effects on offspring's birthweight adjusted for fetal effect", "Maternal effects on offspring's birthweight adjusted for fetal effect", "Maternal effects on offspring's birthweight adjusted for fetal effect", "Maternal effects on offspring's birthweight adjusted for fetal effect",
               "Maternal effects on offspring's birthweight adjusted for fetal effect", "Maternal effects on offspring's birthweight adjusted for fetal effect", "Maternal effects on offspring's birthweight adjusted for fetal effect", "Maternal effects on offspring's birthweight adjusted for fetal effect",
               "Paternal effects on offspring's birthweight", "Paternal effects on offspring's birthweight", "Paternal effects on offspring's birthweight", "Paternal effects on offspring's birthweight",
               "Paternal effects on offspring's birthweight", "Paternal effects on offspring's birthweight", "Paternal effects on offspring's birthweight", "Paternal effects on offspring's birthweight",
               "Paternal effects on offspring's birthweight", "Paternal effects on offspring's birthweight", "Paternal effects on offspring's birthweight", "Paternal effects on offspring's birthweight")


df <- cbind(as.data.frame(betas_), as.data.frame(se_), as.data.frame(pval_), as.data.frame(exposure_), as.data.frame(outcome_), as.data.frame(method_), as.data.frame(MR))

colnames(df) <- c("Beta", "SE", "P", "Exposure", "Outcome", "Method", "MR")

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

tiff("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/PLOTS/2SMR/IS/2SMR_original.tiff", res = 300, height = 2400, width = 6500)
ggplot(df, aes(x=Method, y=as.numeric(Beta), ymin=as.numeric(lower_ci), ymax=as.numeric(upper_ci))) + 
  geom_pointrange(aes(col=Exposure)) +
  geom_errorbar(aes(ymin=as.numeric(lower_ci), ymax=as.numeric(upper_ci),col=Exposure),width=0.5,cex=1)+ 
  geom_hline(aes(fill=Exposure),yintercept =0, linetype=2)+
  xlab('Mendelian Randomization Method')+ ylab("Causal effect with 95% confidence intervals")+
  facet_grid(MR~Outcome, shrink = TRUE, scales = "free") +
  theme(plot.title=element_text(size=12),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12),
        strip.text.y = element_text(hjust=0,vjust = 1))+
   coord_flip()
dev.off()
  
  
  geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = -.9)) +
  scale_fill_manual(values=barCOLS)+
  scale_color_manual(values=dotCOLS)+
  scale_x_discrete(name="Outcome") +
  scale_y_continuous(name = "2SMR effect size with 95% confidence intervals", breaks = c(-0.20, -0.10, 0, 0.10, 0.20), limits = c(-0.25, 0.25)) +
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
