# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: R:percent
#     text_representation:
#       extension: .R
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# %%
rm(list=ls())
libraries = c("dplyr","magrittr","tidyr","reshape2","ggplot2"
              ,"MASS","readr","stats","stringr")
for(x in libraries) { library(x,character.only=TRUE,warn.conflicts=FALSE,quietly=TRUE) }

theme_set(theme_bw())
#setwd("../clade1mpox_SG")

# %%
source("src/reff_model.R")

# %%
source("src/reff_model7.R")

# %%
source("src/reff_model14.R")

# %%
df_SAR_Reff <- read.csv("outputs/outbreakpotential_msm/SAR_m_Reff_4w10d.csv")
head(df_SAR_Reff)

# %%
df_SAR_Reff7 <- read.csv("outputs/outbreakpotential_msm/SAR_m_Reff_4w7d.csv")
head(df_SAR_Reff7)

# %%
df_SAR_Reff14 <- read.csv("outputs/outbreakpotential_msm/SAR_m_Reff_4w14d.csv")
head(df_SAR_Reff14)

# %%
### function of recursive final size equation
final_size_fn = function(sar_ratio,zz){
  N <- 1000000 ### msm size
  RR <- rep(0,zz-1)
  
  for(i in 1:(zz-1)){
    RR[i] <- sar_ratio*df_SAR_Reff$Reff_1[i+1]*(df_SAR_Reff$Infections[i+1]-df_SAR_Reff$Infections[i])
  }
  return(c(sum(RR)/(df_SAR_Reff$Infections[zz])))
}

### solve the equation
final_exposure <- function(sar_ratio){
  x=1500
  X <- rep(0,x)
  for(i in 1:x) {X[i] <- final_size_fn(sar_ratio,i)}   
  
  num <- which(abs(X-1) == min(abs(X-1))) 
  return(df_SAR_Reff$Infections[num])
}

df_final <- mapply(final_exposure, c(1e-5, seq(0.05,10,0.05)))
result_df_final <- cbind(c(1e-5, seq(0.05,10,0.05))*0.1,df_final) %>% as.data.frame()
colnames(result_df_final) <- c("SAR","case")
write.csv(result_df_final, "data/intermediate/SAR_m_final10.csv")


# %%
### function of recursive final size equation
final_size_fn = function(sar_ratio,zz){
  N <- 1000000 ### msm size
  RR <- rep(0,zz-1)
  
  for(i in 1:(zz-1)){
    RR[i] <- sar_ratio*df_SAR_Reff7$Reff_1[i+1]*(df_SAR_Reff7$Infections[i+1]-df_SAR_Reff7$Infections[i])
  }
  return(c(sum(RR)/(df_SAR_Reff7$Infections[zz])))
}

### solve the equation
final_exposure <- function(sar_ratio){
  x=1500
  X <- rep(0,x)
  for(i in 1:x) {X[i] <- final_size_fn(sar_ratio,i)}   
  
  num <- which(abs(X-1) == min(abs(X-1))) 
  return(df_SAR_Reff7$Infections[num])
}

df_final <- mapply(final_exposure, c(1e-5, seq(0.05,10,0.05)))
result_df_final <- cbind(c(1e-5, seq(0.05,10,0.05))*0.1,df_final) %>% as.data.frame()
colnames(result_df_final) <- c("SAR","case")
write.csv(result_df_final, "data/intermediate/SAR_m_final7.csv")

# %%
### function of recursive final size equation
final_size_fn = function(sar_ratio,zz){
  N <- 1000000 ### msm size
  RR <- rep(0,zz-1)
  
  for(i in 1:(zz-1)){
    RR[i] <- sar_ratio*df_SAR_Reff14$Reff_1[i+1]*(df_SAR_Reff14$Infections[i+1]-df_SAR_Reff14$Infections[i])
  }
  return(c(sum(RR)/(df_SAR_Reff14$Infections[zz])))
}

### solve the equation
final_exposure <- function(sar_ratio){
  x=4000
  X <- rep(0,x)
  for(i in 1:x) {X[i] <- final_size_fn(sar_ratio,i)}   
  
  num <- which(abs(X-1) == min(abs(X-1))) 
  return(df_SAR_Reff14$Infections[num])
}

df_final <- mapply(final_exposure, c(1e-5, seq(0.05,10,0.05)))
result_df_final <- cbind(c(1e-5, seq(0.05,10,0.05))*0.1,df_final) %>% as.data.frame()
colnames(result_df_final) <- c("SAR","case")
write.csv(result_df_final, "data/intermediate/SAR_m_final14.csv")

# %%
result_df_final <- read.csv("data/intermediate/SAR_m_final10.csv")
result_df_final7 <- read.csv("data/intermediate/SAR_m_final7.csv")
result_df_final14 <- read.csv("data/intermediate/SAR_m_final14.csv")
result_df_final$R0_10 <- sar_m_reff[1,4] * result_df_final$SAR * 10
result_df_final7$R0_7 <-  sar_m_reff7[1,4] * result_df_final7$SAR * 10
result_df_final14$R0_14 <- sar_m_reff14[1,4] * result_df_final14$SAR * 10
#head(merged_data)

# %%
library(scales)

# %%
# read datasets
sar_m_reff <- read.csv("data/intermediate/SAR_m_Reff_4w10d.csv")
sar_m_final <- read.csv("data/intermediate/SAR_m_final10.csv")
sar_m_reff7 <- read.csv("data/intermediate/SAR_m_Reff_4w7d.csv")
sar_m_final7 <- read.csv("data/intermediate/SAR_m_final7.csv")
sar_m_reff14 <- read.csv("data/intermediate/SAR_m_Reff_4w14d.csv")
sar_m_final14 <- read.csv("data/intermediate/SAR_m_final14.csv")

sar_m_final$case <- as.numeric(sar_m_final$case)
sar_m_final$SAR <- as.numeric(sar_m_final$SAR)
sar_m_final7$case <- as.numeric(sar_m_final7$case)
sar_m_final7$SAR <- as.numeric(sar_m_final7$SAR)
sar_m_final14$case <- as.numeric(sar_m_final14$case)
sar_m_final14$SAR <- as.numeric(sar_m_final14$SAR)
sar_m_reff$SAR <- as.numeric(sar_m_reff$SAR)
sar_m_reff7$SAR <- as.numeric(sar_m_reff7$SAR)
sar_m_reff14$SAR <- as.numeric(sar_m_reff14$SAR)

get_closest_SAR <- function(final_df, reff_df) {
  final_df$closest_SAR <- sapply(final_df$case, function(x) {
    nearest_sar <- reff_df %>%
      filter(abs(Infections - x) == min(abs(Infections - x)))
    
    if (nrow(nearest_sar) == 0) {
      return(NA)
    } else {
      return(nearest_sar$SAR[1])
    }
  })
  return(final_df)
}

# Get the closest SAR
sar_m_final <- get_closest_SAR(sar_m_final, sar_m_reff)
sar_m_final7 <- get_closest_SAR(sar_m_final7, sar_m_reff7)
sar_m_final14 <- get_closest_SAR(sar_m_final14, sar_m_reff14)

label_positions <- seq(min(c(sar_m_final$SAR, sar_m_final7$SAR, sar_m_final14$SAR)), 
                       max(c(sar_m_final$SAR, sar_m_final7$SAR, sar_m_final14$SAR)), 
                       length.out = 5)

# %%
merged_data <- full_join(sar_m_final7, sar_m_final, by = "SAR", suffix = c("_7", "_10")) %>%
  full_join(sar_m_final14, by = "SAR") %>%
  rename(case_14 = case, closest_SAR_14 = closest_SAR)

merged_data <- merged_data %>%
  rename(case_7 = case_7, closest_SAR_7 = closest_SAR_7,
         case_10 = case_10, closest_SAR_10 = closest_SAR_10)

# %%
merged_data <- read.csv("data/intermediate/merged_data.csv")
merged_data$clade1_cR0_7 <- sar_m_reff7[1,4] * merged_data$closest_SAR_7 * 10
merged_data$clade1_cR0_10 <- sar_m_reff[1,4] * merged_data$closest_SAR_10 * 10
merged_data$clade1_cR0_14 <- sar_m_reff14[1,4] * merged_data$closest_SAR_14 * 10
merged_data$clade2_R0_7 <- sar_m_reff7[1,4] * merged_data$SAR * 10
merged_data$clade2_R0_10 <- sar_m_reff[1,4] * merged_data$SAR * 10
merged_data$clade2_R0_14 <- sar_m_reff14[1,4] * merged_data$SAR * 10
head(merged_data)

# %%
write.csv(merged_data,"data/intermediate/merged_data.csv")
