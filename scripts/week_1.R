# EEMB 247 Week 1 Lab
# Numerial Solutions to Ordinary Differential Equations (ODEs)
# date created: April 2, 2026
# date modified: April 2, 2026

#########################################################

install.packages("deSolve")

library(deSolve)
library(tidyverse)
library(here)
library(ggplot2)

##########################################

# basic SIR model

# set parameters
time <- seq(0,5,by=0.01)

beta <- 0.1
gamma <- 1
params <- c(beta, gamma)

s0 <- 99
i0<- 1
r0<- 0
initial <- c(s0,i0,r0)

sir.model <- function(t, x, params){
  s=x[1]
  i=x[2]
  r=x[3]
  beta=params[1]
  gamma=params[2]
  
  dSdt=-beta*s*i
  dIdt=beta*s*i-gamma*i
  dRdt=gamma*i
  
  return(list(c(dSdt,dIdt,dRdt)))
}

results <- lsoda(initial, time, sir.model, params)
colnames(results) <- c("time","s","i","r")


ggplot(results)+
  geom_line(aes(y=s,x=time), color="blue")+
  geom_line(aes(y=r,x=time), color="red")+
  geom_line(aes(y=i,x=time), color="green")+
  theme_bw()

#########################################

# changing initial conditions
time <- seq(0,100,by=0.01)


s0 <- 30
i0<- 1
r0<- 0
initial.adj <- c(s0,i0,r0)

sir.model <- function(t, x, params){
  s=x[1]
  i=x[2]
  r=x[3]
  beta=params[1]
  gamma=params[2]
  
  dSdt=-beta*s*i
  dIdt=beta*s*i-gamma*i
  dRdt=gamma*i
  
  return(list(c(dSdt,dIdt,dRdt)))
}

results <- lsoda(initial.adj, time, sir.model, params)
colnames(results) <- c("time","s","i","r")


ggplot(results)+
  geom_line(aes(y=s,x=time), color="blue")+
  geom_line(aes(y=r,x=time), color="red")+
  geom_line(aes(y=i,x=time), color="green")+
  theme_bw()

############################################

# vaccinate

v <- 0.1

initial.vacc <- c(s0,i0,r0)

time_long <- seq(1,10,by=0.01)

sir.model <- function(t, x, params){
  s=x[1]
  i=x[2]
  r=x[3]
  v=v
  beta=params[1]
  gamma=params[2]
  
  dSdt=-beta*(s-s*v)*i
  dIdt=beta*(s-s*v)*i-gamma*i
  dRdt=gamma*i
  
  return(list(c(dSdt,dIdt,dRdt)))
}

results <- lsoda(initial.vacc, time_long, sir.model, params)
colnames(results) <- c("time","susceptible","infected","recovered")

results.df <- as.data.frame(results)

results.tidy <- results.df |>
  pivot_longer(cols=c("susceptible","infected","recovered"), values_to="humans", names_to ="status")

results.plot <- ggplot(results.tidy, aes(x=time, y=humans, color=status))+
  geom_line()+
  xlim(0,10)+
  ylim(0,50)+
  labs(x="Time (weeks)", y="Human Density", color="Status")+
  scale_color_manual(values=c("blue","red","green"), breaks=c("susceptible", "infected", "recovered"), 
                     labels=c("Susceptible (10% vaccination rate)", "Infected","Recovered"))+
  theme_bw()

results.plot

ggsave(here("figures","week.1.plot.jpg"), results.plot, height=5, width=7, units="in", dpi=300)
