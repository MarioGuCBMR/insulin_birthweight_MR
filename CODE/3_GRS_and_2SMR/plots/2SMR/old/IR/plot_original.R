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

FIadjBMI_fbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/FIadjBMI/FBW/FIadjBMI_FBW_original")
FIadjBMI_mbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/FIadjBMI/MBW/FIadjBMI_MBW_original")
FIadjBMI_fbwadjmbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/FIadjBMI/FBWadjMBW/FIadjBMI_FBWadjMBW_original")
FIadjBMI_mbwadjfbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/FIadjBMI/MBWadjFBW/FIadjBMI_mbwadjfbw_original")
FIadjBMI_father <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/FIadjBMI/Father/FIadjBMI_father_original")

#FI:

FI_fbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/FI/FBW/FI_FBW_original")
FI_mbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/FI/MBW/FI_MBW_original")
FI_fbwadjmbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/FI/FBWadjMBW/FI_FBWadjMBW_original")
FI_mbwadjfbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/FI/MBWadjFBW/FI_mbwadjfbw_original")
FI_father <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/FI/Father/FI_father_original")

#DI:

IR_53_fbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/53_IR_Lotta_SNPs/FBW/IR_53_FBW_original")
IR_53_mbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/53_IR_Lotta_SNPs/MBW/IR_53_MBW_original")
IR_53_fbwadjmbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/53_IR_Lotta_SNPs/FBWadjMBW/IR_53_FBWadjMBW_original")
IR_53_mbwadjfbw <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/53_IR_Lotta_SNPs/MBWadjFBW/IR_53_mbwadjfbw_original")
IR_53_father <- readRDS("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/2SMR/53_IR_Lotta_SNPs/Father/IR_53_father_original")

#################################################
#Let's clean simple mode because nobody likes it#
#################################################

FIadjBMI_fbw <- FIadjBMI_fbw[which(FIadjBMI_fbw$method != "Simple mode"),]
FIadjBMI_fbwadjmbw <- FIadjBMI_fbwadjmbw[which(FIadjBMI_fbwadjmbw$method != "Simple mode"),]
FIadjBMI_mbw <- FIadjBMI_mbw[which(FIadjBMI_mbw$method != "Simple mode"),]
FIadjBMI_mbwadjfbw <- FIadjBMI_mbwadjfbw[which(FIadjBMI_mbwadjfbw$method != "Simple mode"),]
FIadjBMI_father <- FIadjBMI_father[which(FIadjBMI_father$method != "Simple mode"),]

FI_fbw <- FI_fbw[which(FI_fbw$method != "Simple mode"),]
FI_fbwadjmbw <- FI_fbwadjmbw[which(FI_fbwadjmbw$method != "Simple mode"),]
FI_mbw <- FI_mbw[which(FI_mbw$method != "Simple mode"),]
FI_mbwadjfbw <- FI_mbwadjfbw[which(FI_mbwadjfbw$method != "Simple mode"),]
FI_father <- FI_father[which(FI_father$method != "Simple mode"),]

IR_53_fbw <- IR_53_fbw[which(IR_53_fbw$method != "Simple mode"),]
IR_53_fbwadjmbw <- IR_53_fbwadjmbw[which(IR_53_fbwadjmbw$method != "Simple mode"),]
IR_53_mbw <- IR_53_mbw[which(IR_53_mbw$method != "Simple mode"),]
IR_53_mbwadjfbw <- IR_53_mbwadjfbw[which(IR_53_mbwadjfbw$method != "Simple mode"),]
IR_53_father <- IR_53_father[which(IR_53_father$method != "Simple mode"),]

##################################################
#Let's make a dataframe with the betas and the SE#
##################################################

betas_ <- c(FIadjBMI_fbw$b, FI_fbw$b, IR_53_fbw$b, FIadjBMI_fbwadjmbw$b, FI_fbwadjmbw$b, IR_53_fbwadjmbw$b,
           FIadjBMI_mbw$b, FI_mbw$b, IR_53_mbw$b, FIadjBMI_mbwadjfbw$b, FI_mbwadjfbw$b, IR_53_mbwadjfbw$b,
           FIadjBMI_father$b, FI_father$b, IR_53_father$b)

se_ <- c(FIadjBMI_fbw$se, FI_fbw$se, IR_53_fbw$se, FIadjBMI_fbwadjmbw$se, FI_fbwadjmbw$se, IR_53_fbwadjmbw$se,
           FIadjBMI_mbw$se, FI_mbw$se, IR_53_mbw$se, FIadjBMI_mbwadjfbw$se, FI_mbwadjfbw$se, IR_53_mbwadjfbw$se,
           FIadjBMI_father$se, FI_father$se, IR_53_father$se)

pval_ <-  c(FIadjBMI_fbw$pval, FI_fbw$pval, IR_53_fbw$pval, FIadjBMI_fbwadjmbw$pval, FI_fbwadjmbw$pval, IR_53_fbwadjmbw$pval,
            FIadjBMI_mbw$pval, FI_mbw$pval, IR_53_mbw$pval, FIadjBMI_mbwadjfbw$pval, FI_mbwadjfbw$pval, IR_53_mbwadjfbw$pval,
            FIadjBMI_father$pval, FI_father$pval, IR_53_father$pval)

exposure_ <-  c("FIadjBMI", "FIadjBMI", "FIadjBMI", "FIadjBMI", "FI", "FI", "FI", "FI", "IR_53", "IR_53", "IR_53", "IR_53",
                "FIadjBMI", "FIadjBMI", "FIadjBMI", "FIadjBMI", "FI", "FI", "FI", "FI", "IR_53", "IR_53", "IR_53", "IR_53",
                "FIadjBMI", "FIadjBMI", "FIadjBMI", "FIadjBMI", "FI", "FI", "FI", "FI", "IR_53", "IR_53", "IR_53", "IR_53",
                "FIadjBMI", "FIadjBMI", "FIadjBMI", "FIadjBMI", "FI", "FI", "FI", "FI", "IR_53", "IR_53", "IR_53", "IR_53",
                "FIadjBMI", "FIadjBMI", "FIadjBMI", "FIadjBMI", "FI", "FI", "FI", "FI", "IR_53", "IR_53", "IR_53", "IR_53")

method_ <-  c("Egger-FIadjBMI", "Weigthed_Median_FIadjBMI", "IVW_FIadjBMI", "Weighted_Mode_FIadjBMI", "Egger-FI", "Weighted_Median_FI", "Egger_FI", "Weighted_Mode_FI", "Egger_DI", "Weighted_Median_DI", "IVW_DI", "Weighted_Mode_DI",
              "Egger-FIadjBMI", "Weigthed_Median_FIadjBMI", "IVW_FIadjBMI", "Weighted_Mode_FIadjBMI", "Egger-FI", "Weighted_Median_FI", "Egger_FI", "Weighted_Mode_FI", "Egger_DI", "Weighted_Median_DI", "IVW_DI", "Weighted_Mode_DI",
              "Egger-FIadjBMI", "Weigthed_Median_FIadjBMI", "IVW_FIadjBMI", "Weighted_Mode_FIadjBMI", "Egger-FI", "Weighted_Median_FI", "Egger_FI", "Weighted_Mode_FI", "Egger_DI", "Weighted_Median_DI", "IVW_DI", "Weighted_Mode_DI",
              "Egger-FIadjBMI", "Weigthed_Median_FIadjBMI", "IVW_FIadjBMI", "Weighted_Mode_FIadjBMI", "Egger-FI", "Weighted_Median_FI", "Egger_FI", "Weighted_Mode_FI", "Egger_DI", "Weighted_Median_DI", "IVW_DI", "Weighted_Mode_DI",
              "Egger-FIadjBMI", "Weigthed_Median_FIadjBMI", "IVW_FIadjBMI", "Weighted_Mode_FIadjBMI", "Egger-FI", "Weighted_Median_FI", "Egger_FI", "Weighted_Mode_FI", "Egger_DI", "Weighted_Median_DI", "IVW_DI", "Weighted_Mode_DI")

MR <-  c(FIadjBMI_fbw$method, FIadjBMI_fbw$method, FIadjBMI_fbw$method, FIadjBMI_fbw$method, FIadjBMI_fbw$method,
         FIadjBMI_fbw$method, FIadjBMI_fbw$method, FIadjBMI_fbw$method, FIadjBMI_fbw$method, FIadjBMI_fbw$method,
         FIadjBMI_fbw$method, FIadjBMI_fbw$method, FIadjBMI_fbw$method, FIadjBMI_fbw$method, FIadjBMI_fbw$method)


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

tiff("C:/Users/zlc436/Desktop/2SMR_4_Hermina/OUTPUT/PLOTS/2SMR/IR/2SMR_original.tiff", res = 300, height = 2400, width = 6500)
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
