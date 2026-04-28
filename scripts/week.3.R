# EEMB 247
# Eva Juengling Bean
# Week 3: Stability Analyses
# created 4/21/2026
# last modified: 4/28/2026
############################################

rm(list=ls())
# install.packages("rootSolve")
############################################

#load libraries

library(tidyverse)
library(here)
library(janitor)
library(deSolve)
library(rootSolve)
library(usethis)
library(ggplot2)

# load in data

frog.data <- read.csv(here("data","frog_data.csv"))

# simulation equilibrium

IN=10
dN=0.2

params <- c(IN,dN)

prey.model <- function(N, state, params){
  N=state[1]
  IN = params[1]
  dN = params[2]
  
  dNdt=IN-dN*N
  
  return(list(dNdt))
}

N0=50
times <- seq(0,100,by=0.1)

results <- lsoda(N0, times, prey.model, params)
colnames(results) <- c("time", "N")

results.df <- data.frame(results)

results.plot <- ggplot(results.df, aes(x=time, y=N))+
  geom_line()+
  theme_bw()

results.plot

# numerical solution

rhs.prey <- function(N, IN, dN){
  IN-dN*N
}

prey.equil <- uniroot.all(rhs.prey, c(0,10), IN=IN, dN=dN)

prey.equil

IN.all <- seq(0,20, by=0.1)

prey.equilibrium <- rep(0, length(IN.all))

for (i in 1:length(IN.all)){

  prey.equilibrium[i] = uniroot.all(rhs.prey, c(0,200), IN=IN.all[i], dN=dN)
}

equil.df <- data.frame(IN.all, prey.equilibrium)

equil.plot <- ggplot(equil.df)+
  geom_line(aes(x=IN.all, y=prey.equilibrium))+
  theme_bw()

equil.plot

########################################################

# Predator and prey
# simulation solution

IN <- 10
dN <- 0.2
a <- 0.01
c <- 0.5
dP <- 0.1

N0 <- 2
P0 <- 2

params <- c(IN, dN, a, c, dP)
initial.states <- c(N0, P0)

predator.prey <- function(t, states, params){
  IN = params[1]
  dN = params[2]
  a = params[3]
  c = params[4]
  dP = params[5]
  
  N = states[1]
  P = states[2]
  
  dNdt = IN-dN*N-a*N*P
  dPdt = c*a*N*P-dP*P
  
  return(list(c(dNdt,dPdt)))
}

times <- seq(0,100, by=0.1)

results <- lsoda(initial.states, times, predator.prey, params)
colnames(results)<- c("time", "prey", "predator")

results.data <- data.frame(results)

plot <- ggplot(results.data)+
  geom_line(aes(x=time, y=prey), color="darkred")+
  geom_line(aes(x=time, y=predator), color="forestgreen")+
  theme_bw()

plot

parameters <- c(IN=10, dN=0.2,a=0.01, c=0.5, dP=0.1)
x <- c(N=10, P=10)

predator.prey.root <- function(x, parameters){
  IN = params[1]
  dN = params[2]
  a = params[3]
  c = params[4]
  dP = params[5]
  
  N = x[1]
  P = x[2]
  
  F1 = IN-dN*N-a*N*P
  F2 = c*a*N*P-dP*P
  
  return(c(F1, F2))
  }

eqm <- multiroot(predator.prey.root, x, parms = parameters)

eqm$root

###############################################################

# Nicholson-Bailey host-parasitoid model
# simulation

max.time <- 50
initial.p <- 11
initial.h <- 27

R <- 2
a <- 0.05
c <- 0.9

p <- rep(0, max.time+1)
h <- rep(0, max.time+1)
p[1] <- initial.p
h[1] <- initial.h

for (t in 1:max.time){
  h[t+1] = R*h[t]*exp(-a*p[t])
  p[t+1]=c*h[t]*(1-exp(-a*p[t]))
}

times <- seq(0,max.time)

df <- data.frame(p, h, times)

plot <- ggplot(df)+
  geom_line(aes(x=times, y=p), color="forestgreen")+
  geom_line(aes(x=times, y=h), color="darkred")+
  theme_bw()

plot # plot has big spikes and dies out--> no stable equilibrium

################################################################

# numerical solution

parameters <- c(R=R, a=a, c=c)

n.b.eqm <- function(x, parms){
  h=x[1]
  p=x[2]
  
  R=parms[1]
  a=parms[2]
  c=parms[3]
  
  F1=h-R*h*exp(-a*p)
  F2=p-c*h*(1-exp(-a*p))

  return(c(F1,F2))
}

x <- c(h=10, p=10)

eqm.solve <- multiroot(n.b.eqm,x,parms=parameters)

eqm.solve$root # solutions the same as calculated!
