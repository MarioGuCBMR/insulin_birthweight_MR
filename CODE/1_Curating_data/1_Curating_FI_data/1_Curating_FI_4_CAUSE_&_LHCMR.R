##############
#INTRODUCTION#
##############

#This code is to clean the Lagou et al 2021 FI data:

###################
#Loading Libraries#
###################

library(data.table)
library(tidyverse)

##############
#INTRODUCTION#
##############

#This is a code to solve the build issues with FI data:

#Probably this won't affect the results.
#But better be save than sorry.

###########
#libraries#
###########

library(data.table)
library(tidyverse)

###########
#Functions#
###########

esv_parser <- function(rsid){
  
  #This function removes the weirds ens SNPs...
  
  tmp <- unlist(strsplit(rsid, "esv"))
  
  if(length(tmp)<2){
    
    return(rsid)
    
  }
  
}

hg_parser <- function(rsid){
  
  #This function removes the weirds hg SNPs...
  
  tmp <- as.character(unlist(strsplit(rsid, "hg18")))
  
  if(length(tmp)==2){
    
    tmp_  <- as.character(unlist(strsplit(tmp, "_")))
    
    chr_pos <- paste("chr", tmp_[2], sep = "")
    
    chr_pos_ <- paste(chr_pos, tmp_[3], sep = ":")
    
    return(chr_pos_)
    
  } else {
    
    return(rsid)
    
  }
  
}

chr_parser <- function(chr_pos){
  #This function gets the chromosome of a variant in chr: format.
  
  check_separator <- str_detect(chr_pos, "_")
  
  if(check_separator == FALSE){
    
    sep = ":"
    
  } else {
    
    sep = "_"
    
  }
  
  tmp <- strsplit(chr_pos, sep)[[1]][1]
  
  return(tmp)
  
}

pos_parser <- function(chr_pos){
  #This function generates the position for liftover format.
  
  check_separator <- str_detect(chr_pos, "_")
  
  if(check_separator == FALSE){
   
    sep = ":"
     
  } else {
    
    sep = "_"
    
  }
  
  tmp <- strsplit(chr_pos, sep)[[1]][2]
  
  tmp_1 <- as.numeric(as.character(tmp)) - 1
  tmp_2 <- as.numeric(as.character(tmp)) + 1
  
  tmp_end <- paste(tmp_1, tmp_2, sep = "-")
  
  return(tmp_end)
  
}

weird_allele_parser <- function(chr_pos){
  #Cleans the weird alleles which have their chromosomes in their position too.
  
  #1. first we split the chr_pos in chr and pos:
  
  chr_ <- strsplit(chr_pos, ":")[[1]][1]
  pos_ <- strsplit(chr_pos, ":")[[1]][2]
  
  #2. let's check whether we need to prune the string 1 or 2 characters depending on the chromosome.
  
  check = nchar(chr_) #5 if it is >= 10. 4 < 10
  
  if(check == 5){
    
    pos_curated <- str_sub(pos_, 3)
    
    chr_pos_curated <- paste(chr_, ":", pos_curated, sep = "")
    
    return(chr_pos_curated)
    
  } else { #there is only one other option so...
    
    pos_curated <- str_sub(pos_, 2)
    
    chr_pos_curated <- paste(chr_, ":", pos_curated, sep = "")
    
    return(chr_pos_curated)
    
  }

}
  


liftover_parser <- function(chr_pos){
  
  chr <-  strsplit(chr_pos, ":")[[1]][1]
  
  pos <- strsplit(chr_pos, "-")[[1]][2]
  
  pos_corrected <- as.numeric(as.character(pos)) - 1
  
  chr_pos_37 <- paste(chr, as.character(pos_corrected), sep = ":")
  
  return(chr_pos_37)
  
}

parse_rsid <- function(rsid){
  
  rsid_ <- strsplit(rsid, ":")[[1]][1] #this works for those that do not present the weird formatting :)
  
  return(rsid_)
  
}

#################
#Loading FI data#
#################

#This FI data is the one that we have from Hermina's project.
#This data has rsid that are mostly from build 37. Some are from build 36, but since we have no chromosome and position
#there is no way for me to retrieve them.

#1) There is not key line getter from build 36 to 37.

#2) Phenoscanner only erases those that are not in the new Build. Hence, they would all be removed either way.
#Since they are all removed..., then it does not matter if they are there or not. They would just not match with the new build.

#3) Regarding the merging of RSIDs, we could control that with build 37 to 38, but not for 36 to 37.
#Using phenoscanner the only thing that we are taking into account are those in 36 that are the same as in 37
#but they are merged in 38.

#4) Hence, with RSID that is as much as we can do. This does not affect the preliminary results: Lotta results with
#RSID from build 36 are exactly the same as with Lotta's.
#Hence, I think we are gonna be okay.

#5) Nonetheless, that does not mean that we cannot work with the variants that are weird. We have some variants that present chromosome and position.
#These variants can be changed to RSID and we would be missing SNPs if we are forgetting them.
#Hence, it is essential to run these.
#For Hermina's project we already did it. So, let's go.

FI_original <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/RAW_DATA/Glycemic_traits_raw_data/FI_combined_1000G_density.txt.gz")

######
#MAIN#
######

#We know this data by heart and we know that we have to include only those SNPs that are GWAS_inserted.
#The others do not present data from European sources.

FI_original <- FI_original[which(FI_original$source != "SSIMP"),]

########
#STEP 1#
########

#First we are going to remove the low sample sizes:

#We decided to go for those with a sample size below 10.000 N

effect_sample_fi <- 10000

FI_original <- FI_original[which(FI_original$n > effect_sample_fi),]

summary(FI_original$r2.pred) #the info here is 1 since they were not imputed. 
summary(FI_original$eaf) #Some of them need to be removed:

FI_original <- FI_original[which(FI_original$eaf > 0.01),]

summary(FI_original$eaf) #PERFECT.

FI_original <- FI_original[order(FI_original$rsid),]

head(FI_original, 1000) #we do have some of them in chr_pos format only. The readme does not say in which build.

#Let's check what kind of data do we have:

FI_original_chr <- FI_original[which(str_detect(FI_original$rsid, "rs") == FALSE),] #9860 variants have chr_pos
FI_original_rsid <- FI_original[which(str_detect(FI_original$rsid, "rs") == TRUE),] #here we have the rest..., but do we only have rsid???

FI_original_rsid <- FI_original_rsid[order(FI_original_rsid$rsid),]

head(FI_original_rsid, 1000) #yess, all is good.
head(FI_original_chr, 1000) #only chr_pos here
tail(FI_original_chr, 100) #we have some variants with hg format!

FI_original_chr_clean <-  FI_original_chr[which(str_detect(FI_original_chr$rsid, "chr") == TRUE),] #9698 have chr.

#What about the rest?

FI_original_not_chr_or_rsid <-  FI_original_chr[which(str_detect(FI_original_chr$rsid, "chr") == FALSE),] #162 have chr.

#Let's check what do we have here:

FI_original_not_chr_or_rsid <- FI_original_not_chr_or_rsid[order(FI_original_not_chr_or_rsid$rsid),]

head(FI_original_not_chr_or_rsid)
tail(FI_original_not_chr_or_rsid)

#################################################
#CHECKING CHR POS ONES: IN WHICH BUILD ARE THEY?#
#################################################

#We can know this with phenoscanner. PhenoScaner only works with builds 37 and 38.
#Meaning that if most are not found..., they are most probably from build hg18 too.

random_variants <- sample(FI_original_chr_clean$rsid, 100) #Let's check with random variants.

ps_results <- phenoscanner::phenoscanner(random_variants)

ps_results_end <- ps_results$snps #only one variant was correctly converted: chr2:203739475.

FI_original[which(FI_original$rsid == "chr2:203739475"),] #A/T in FI data.
ps_results_end #A/G in phenoscanner.

#Let's check in dbsnp if it is a triallelic SNP... IT IS NOT.
#It is a bi-allelic variant.
#Indicating that it was just a random match.
#This variant IS NOT in build 19.

#Let's do the conversion from 18 to 19 and check whether the alleles match then.

bed_18_to_19 <- as.data.frame(read.table("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/FI/test_18_to_19.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

#chr2:204031229-204031231

#Let's check this in dbsnp: 2:204031230
#The SNP is bi-allelic! Now the alleles match.
#I think this is hard proof: all the variants are in build hg18.

##############################################
#PREPARING DATA TO BE TRANSFORMED IN LIFTOVER#
##############################################

#First we need to get the data for liftover ready from the different sources.
#We have the data in hg format and the chr_clean data.

#Let's clean them individually:

head(FI_original_not_chr_or_rsid) #for this ones it is really easy:

#In these conversions we need to be wary how R treats the numbers so...

options(scipen = 100) #important to not make liftover crazy!

FI_original_not_chr_or_rsid$chr_pos_18 <- as.character(unlist(sapply(FI_original_not_chr_or_rsid$rsid, hg_parser)))
FI_original_not_chr_or_rsid$pos_18_4_liftover <- as.character(unlist(sapply(FI_original_not_chr_or_rsid$chr_pos_18, pos_parser)))
FI_original_not_chr_or_rsid$chr_4_liftover <- as.character(unlist(sapply(FI_original_not_chr_or_rsid$chr_pos_18, chr_parser)))
FI_original_not_chr_or_rsid$chr_pos_18_4_liftover <- paste(FI_original_not_chr_or_rsid$chr_4_liftover, FI_original_not_chr_or_rsid$pos_18_4_liftover, sep = ":")

FI_original_not_chr_or_rsid <- FI_original_not_chr_or_rsid %>%
  select(-c(chr_4_liftover, pos_18_4_liftover))

FI_original_not_chr_or_rsid <- FI_original_not_chr_or_rsid[order(FI_original_not_chr_or_rsid$rsid),]

head(FI_original_not_chr_or_rsid) #everything seems to be in order.
tail(FI_original_not_chr_or_rsid) #everything seems to be in order.

#Let's check the ranges just in case:

FI_original_not_chr_or_rsid <- FI_original_not_chr_or_rsid[order(FI_original_not_chr_or_rsid$chr_pos_18_4_liftover),]

View(FI_original_not_chr_or_rsid$chr_pos_18_4_liftover)

#Now the same for the clean ones, though in this case we just need to use pos_parser to get the data ready for liftover:

FI_original_chr_clean$chr_pos_18 <- FI_original_chr_clean$rsid
FI_original_chr_clean$pos_18_4_liftover <- as.character(unlist(sapply(FI_original_chr_clean$chr_pos_18, pos_parser)))
FI_original_chr_clean$chr_4_liftover <- as.character(unlist(sapply(FI_original_chr_clean$chr_pos_18, chr_parser)))
FI_original_chr_clean$chr_pos_18_4_liftover <- paste(FI_original_chr_clean$chr_4_liftover, FI_original_chr_clean$pos_18_4_liftover, sep = ":")

FI_original_chr_clean <- FI_original_chr_clean %>%
  select(-c(chr_4_liftover, pos_18_4_liftover))

FI_original_chr_clean <- FI_original_chr_clean[order(FI_original_chr_clean$rsid),]
 
head(FI_original_chr_clean) #everything seems to be in order.
tail(FI_original_chr_clean) #everything seems to be in order.

FI_original_4_liftover <- rbind(FI_original_chr_clean, FI_original_not_chr_or_rsid)

#NOTE: some conversions have failed and we do not know why since the numbers reported in liftover do not make sense.
#I will try to investigate that by separting the liftover in 1000:

fwrite(as.data.frame(FI_original_4_liftover$chr_pos_18_4_liftover), "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/FI/FI_build38.txt", sep = " ", col.names = FALSE)

#29 conversions failed..., nothing that we can do there!
#We are going to split the data according to the failed and not failed and go again:

failures <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/FI/Failed_Conversions.txt")

#We have very weird failures and that is because some rsid cannot be properly converted.
#I have the theory that these indicate that the position is wrong and that is larger than the chromosome for some reason.
#Hence what I am going to do is check which variants have positions larger than their chromosome size.
#There shouldn't be any..., but I am sure they exist.

FI_original_4_liftover$pos_18_4_liftover <- as.character(unlist(sapply(FI_original_4_liftover$chr_pos_18, pos_parser)))
FI_original_4_liftover$chr_4_liftover <- as.character(unlist(sapply(FI_original_4_liftover$chr_pos_18, chr_parser)))

#Just by looking at the first example I already know what is going on. 
#These big position actually have their chromosome as their first two digits!!

#To properly identify with these variants we are going to check in which chromosomes we have the difficult SNPs, and identify them by chromosome size.

size_1 <- 	248956422
size_2 <- 	242193529
size_3 <- 	198295559
size_4 <- 	190214555
size_5 <- 181538259
size_6 <- 	170805979 
size_7 <- 159345973
size_8 <- 	145138636
size_9 <- 	138394717
size_10 <- 133797422 
size_11 <- 135086622 
size_12 <- 133275309
size_13 <- 114364328
size_14 <- 107043718
size_15 <- 101991189
size_16 <- 90338345
size_17 <- 83257441
size_18 <- 80373285
size_19 <- 58617616
size_20 <- 64444167
size_21 <- 46709983
size_22 <- 50818468

size_vect <- c(size_1, size_2, size_3, size_4, size_5, size_6, size_7, size_8, size_9, size_10,
               size_11, size_12, size_13, size_14, size_15, size_16, size_17, size_18, size_19, size_20,
               size_21, size_22)


#Let' save the weird SNPs and use them to clean the data and liftover it.

weird_snps <- c()

for(index_ in seq(1, length(FI_original_4_liftover$rsid))){
  
  chr_pos = FI_original_4_liftover$chr_pos_18[index_]
  chr = FI_original_4_liftover$chr_4_liftover[index_]
  pos = FI_original_4_liftover$pos_18_4_liftover[index_]
  
  for(chromosome in seq(1,22)){
    
    chr_ <- paste("chr", as.character(chromosome), sep = "")
    
    #print(chr_)
    
    if(chr == chr_){
      
      final_size <- sum(size_vect[1:chromosome])
      
      #print(chr)
      
      pos_ <- strsplit(pos, "-")[[1]][1]
      
      if(as.numeric(pos_) > as.numeric(final_size)){
        
        #print(pos_)
        #print(final_size)
        weird_snps <- c(weird_snps, chr_pos)
        
        #return(chr_pos)
        
      }
      
    }
    
  }
  
}



weird_snps <- as.character(unlist(sapply(FI_original_4_liftover$chr_pos_18, weird_allele_parser, chr = FI_original_4_liftover$chr_4_liftover, pos = FI_original_4_liftover$pos_18_4_liftover)))

#IT IS JUST AS I EXPECTED: THEIR FIRST DIGIT IS ACTUALLY THEIR CHROMOSOME.

#Let's divide the dataframe according to these variants:

FI_original_4_liftover_OK <- FI_original_4_liftover[which(!(FI_original_4_liftover$chr_pos_18%in%weird_snps)),]
FI_original_4_liftover_weird <- FI_original_4_liftover[which(FI_original_4_liftover$chr_pos_18%in%weird_snps),]

#We need to clean the chr_pos_18 and rsid variants first, then run the rest

FI_original_4_liftover_weird$rsid <- as.character(unlist(sapply(FI_original_4_liftover_weird$rsid, weird_allele_parser)))
FI_original_4_liftover_weird$chr_pos_18 <- as.character(unlist(sapply(FI_original_4_liftover_weird$chr_pos_18, weird_allele_parser)))
FI_original_4_liftover_weird$pos_18_4_liftover <- as.character(unlist(sapply(FI_original_4_liftover_weird$chr_pos_18, pos_parser)))
FI_original_4_liftover_weird$chr_4_liftover <- as.character(unlist(sapply(FI_original_4_liftover_weird$chr_pos_18, chr_parser)))
FI_original_4_liftover_weird$chr_pos_18_4_liftover <- paste(FI_original_4_liftover_weird$chr_4_liftover, FI_original_4_liftover_weird$pos_18_4_liftover, sep = ":")

#Some of the positions start with 0s and I have no clue if that makes sense.
#My code naturally removes the 0s, as expected.
#Maybe this won't make any issue... to test for this I am gonna run chr10:12635324-12635326 on liftover.
#If it does not make sense with dbsnp... I am going to drop these variants.
#They are just too confusing.

test_weird <- as.data.frame(read.table("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/FI/weird_snp_test.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

#This didn't work out as I expected we are run liftover with the usual data:

fwrite(as.data.frame(FI_original_4_liftover_OK$chr_pos_18_4_liftover), "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/FI/FI_build38_no_weird.txt", sep = " ", col.names = FALSE)

#One of them is deleted, which is probably in the failures file from above.
#I will add clean the data:

weird_snps <- c("chr2:21043776", weird_snps)

FI_original_4_liftover_OK <- FI_original_4_liftover[which(!(FI_original_4_liftover$chr_pos_18%in%weird_snps)),]
FI_original_4_liftover_weird <- FI_original_4_liftover[which(FI_original_4_liftover$chr_pos_18%in%weird_snps),]

fwrite(as.data.frame(FI_original_4_liftover_OK$chr_pos_18_4_liftover), "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/FI/FI_build38_no_weird.txt", sep = " ", col.names = FALSE)

######################################
#Let's recover the data from liftover#
######################################

FI_liftover_results <- read.table("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/FI/FI_build37_raw.bed", quote = "", sep = "\t")

FI_original_4_liftover_OK$chr_pos_37 <- as.character(unlist(sapply(FI_liftover_results$V1, liftover_parser)))

#And now let's merge it with those that failed:

FI_original_4_liftover_weird$chr_pos_37 <- "-"

FI_original_4_liftover_recovered <- rbind(FI_original_4_liftover_OK, FI_original_4_liftover_weird)

##########################################################################
#Next step: let's recover the chr_pos of the variants that only have rsID#
##########################################################################

#The first step is to obtain as many matches as possible by using BMI data from Pulit et al:

bmi <- fread("C:/Users/zlc436/Desktop/PULIT_data_reservoir/Combined/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz")

#We need to be careful, some SNPs from Pulit et al do not present correct chr_pos.
#The mistake is in those without INFO data, thus, the ones in GIANT data.
#We are gonna remove those.

summary(bmi$INFO) 
summary(bmi$Freq_Tested_Allele) 
summary(bmi$CHR) 
summary(bmi$N) 

bmi <- bmi[which(is.na(bmi$INFO) == FALSE),]

summary(bmi$INFO) #all good
summary(bmi$Freq_Tested_Allele) #all good
summary(bmi$CHR) #all good
summary(bmi$N) #all good

#It seems all is OK.
#Now we need to clean the RSIDs, because the have the effect allele and the other allele.
#Which is actually fine.

bmi$RSID <- as.character(unlist(sapply(bmi$SNP, parse_rsid)))

#Let's match the data:

bmi_match <- bmi[which(bmi$RSID%in%FI_original_rsid$rsid),]
fi_match <- FI_original_rsid[which(FI_original_rsid$rsid%in%bmi$RSID),]

#We have some duplicates in BMI, but we do not care about those!
#We only want the positions, the alleles are not the focus here:

bmi_match_no_dupl <- bmi_match[which(duplicated(bmi_match$RSID) == FALSE),] #PERFECT MATCH WITH FI_match.

fi_match <- fi_match[order(match(fi_match$rsid, bmi_match_no_dupl$RSID)),]

length(which(fi_match$rsid == bmi_match_no_dupl$RSID)) #all of them. AWESOME.

fi_match$chr_pos_37 <- paste("chr", bmi_match_no_dupl$CHR, ":", bmi_match_no_dupl$POS, sep="")

############################################################
#The rest are gonna be recovered in SNPnexus + phenoscanner#
############################################################

fi_mismatch <- FI_original_rsid[which(!(FI_original_rsid$rsid%in%bmi$RSID)),] #27.000 variants.

fi_4_nexus <- fi_mismatch

fi_4_nexus$dbsnp <- "dbsnp"

fi_4_nexus_clean <- fi_4_nexus %>%
  select(dbsnp, rsid)

#That is perfect. 

fwrite(fi_4_nexus_clean, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/FI/FI_combined_4_nexus.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

#####################
#Let's retreive them#
#####################

fi_nexus_results <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/FI/FI_Nexus_results.txt")

fi_mismatch_solved <- fi_mismatch[which(fi_mismatch$rsid%in%fi_nexus_results$dbSNP),]
fi_nexus_results <- fi_nexus_results[order(match(fi_nexus_results$dbSNP, fi_mismatch_solved$rsid)),]

which(fi_nexus_results$dbSNP != fi_mismatch_solved$rsid) #they all match.

fi_nexus_results$chr_ <-  paste("chr", fi_nexus_results$Chromosome, sep = "")
fi_nexus_results$chr_pos <- paste(fi_nexus_results$chr_, fi_nexus_results$Position, sep = ":")

#worked perfectly.

fi_mismatch_solved$chr_pos_37 <- fi_nexus_results$chr_pos

#########
#PERFECT#
#########

#Only some 3000 variants to go:

fi_mismatch_failed <- fi_mismatch[which(!(fi_mismatch$rsid%in%fi_nexus_results$dbSNP)),]

#These are quite many, we can query up to 10.000 SNPs in an hour... so it can be done without an issue.
#A bit of a pain, but we will have to do line by line since, if not, the API gets crazy and shuts down because we are asking
#for a bit too much. Sorry my fellow bioinformaticians, you are gonna see some bad coding here.

ps_results_1 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[1:100])$snp
ps_results_2 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[101:200])$snp
ps_results_3 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[201:300])$snp
ps_results_4 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[301:400])$snp
ps_results_5 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[401:500])$snp
ps_results_6 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[501:600])$snp
ps_results_7 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[601:700])$snp
ps_results_8 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[701:800])$snp
ps_results_9 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[801:900])$snp
ps_results_10 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[901:1000])$snp
ps_results_11 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[1001:1100])$snp
ps_results_12 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[1101:1200])$snp
ps_results_13 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[1201:1300])$snp
ps_results_14 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[1301:1400])$snp
ps_results_15 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[1401:1500])$snp
ps_results_16 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[1501:1600])$snp
ps_results_17 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[1601:1700])$snp
ps_results_18 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[1701:1800])$snp
ps_results_19 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[1801:1900])$snp
ps_results_20 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[1901:2000])$snp
ps_results_21 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[2001:2100])$snp
ps_results_22 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[2101:2200])$snp
ps_results_23 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[2201:2300])$snp
ps_results_24 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[2301:2400])$snp
ps_results_25 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[2401:2500])$snp
ps_results_26 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[2501:2600])$snp
ps_results_27 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[2601:2700])$snp
ps_results_28 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[2701:2800])$snp
ps_results_29 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[2801:2900])$snp
ps_results_30 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[2901:3000])$snp
ps_results_31 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[3001:3100])$snp
ps_results_32 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[3101:3200])$snp
ps_results_33 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[3201:3300])$snp

ps_results_34 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[3301:3400])$snp
ps_results_35 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[3401:3500])$snp
ps_results_36 <- phenoscanner::phenoscanner(fi_mismatch_failed$rsid[3501:3509])$snp

ps_results <- rbind(ps_results_1, ps_results_2, ps_results_3, ps_results_4, ps_results_5,
                    ps_results_6, ps_results_7, ps_results_8, ps_results_9, ps_results_10,
                    ps_results_11, ps_results_12, ps_results_13, ps_results_14, ps_results_15,
                    ps_results_16, ps_results_17, ps_results_18, ps_results_19, ps_results_20,
                    ps_results_21, ps_results_22, ps_results_23, ps_results_24, ps_results_25,
                    ps_results_26, ps_results_27, ps_results_28, ps_results_29, ps_results_30,
                    ps_results_31, ps_results_32, ps_results_33, ps_results_34, ps_results_35, ps_results_36)
                    

###################################################
#I am getting pretty sick of this, ain't gonna lie#
###################################################

fi_mismatch_failed_solved <- fi_mismatch_failed[which(fi_mismatch_failed$rsid%in%ps_results$snp),]

ps_results <- ps_results[order(match(ps_results$snp, fi_mismatch_failed_solved$rsid)),]

length(which(ps_results$snp != fi_mismatch_failed_solved$rsid))
length(which(ps_results$snp == fi_mismatch_failed_solved$rsid))

fi_mismatch_failed_solved$chr_pos_37 <- ps_results$hg19_coordinates

#The rest of SNPs are going to just be cleaned left without chr_pos.
#We did recover most of them, so great!

fi_mismatch_failed_unsolved <- fi_mismatch_failed[which(!(fi_mismatch_failed$rsid%in%ps_results$snp)),]

fi_mismatch_failed_unsolved$chr_pos_37 <- "-"

##################################
#Now we merge all the data needed#
##################################

#First let's merge the data from PhenoScanner:

fi_check_on_ps <- rbind(fi_mismatch_failed_solved, fi_mismatch_failed_unsolved)

#Now we merge this one with those recovered with dbsnp:

fi_check_on_dbsnp <- rbind(fi_check_on_ps, fi_mismatch_solved)

#Now we merge this one with those recovered with Pulit data:

fi_check_on_pulit <- rbind(fi_check_on_dbsnp, fi_match)

#Now we are going to clean the data on those that we liftovered to be able to merge the data:

head(FI_original_4_liftover_recovered)
head(fi_check_on_pulit)

fi_check_on_liftover <- FI_original_4_liftover_recovered %>%
  dplyr::select(-c("chr_pos_18_4_liftover", "chr_pos_18", "pos_18_4_liftover", "chr_4_liftover"))

#And now we merge them:

fi_curated <- rbind(fi_check_on_pulit, fi_check_on_liftover) #same amount of rows as the original. WE GOT EM.

#We got them boys!!

fwrite(fi_curated, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/FI/FI_combined_curated.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

check <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/FI/FI_combined_curated.txt")

##################################################################################################################
#FINAL NOTE: we will remove the MHC by overlapping the curated lipid data which REALLY has the MHC region removed#
##################################################################################################################