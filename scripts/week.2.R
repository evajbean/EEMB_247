# EEMB 247
# Eva Juengling Bean
# Week 2 -- Discrete and Logistic Time
#########################################

rm(list=ls())

library(deSolve)
library(tidyverse)
library(janitor)
library(here)
library(ggplot2)

# set parameters

a11 <- 1
a12 <- 0.5
a21 <- 0.5
a22 <- 1

A <- matrix(c(a11,a12,a21,a22),nrow=2, ncol=2, byrow=TRUE)
A

eigen(A)

# numerical solution

initial.x <- 0
initial.y <- 1
max.time <- 10

x <- rep(0,max.time+1)
y <- rep(0, max.time + 1)
 x[1]= initial.x
 y[1]=initial.y
 
 for (t in 1:max.time){
   x[t+1]=a11*x[t]+a12*y[t]
   y[t+1]=a21*x[t]+a22*y[t]
 }
 
times <- seq(0,max.time)

df <- data.frame(x, y, times)

discrete.plot <- ggplot(df)+
  geom_line(aes(x=times, y=x), color="blue")+
  geom_line(aes(x=times, y=y), color="red")+
  xlab("Time")+
  ylab("X or Y")+
  theme_bw()

discrete.plot

# analytical solution:
# Constants = [0.5
#              0.5]
x_analytical <-  0.5*(1.5^times)-0.5*(0.5^times)
y_analytical <- 0.5*(1.5^times)+0.5*(0.5^times)

df <- df |>
  mutate("x.analytical"=x_analytical, "y.analytical"=y_analytical)

discrete.plot <- ggplot(df)+
  geom_line(aes(x=times, y=x), color="blue")+
  geom_line(aes(x=times, y=y), color="red")+
  geom_point(aes(x=times, y=x.analytical, color="skyblue"))+
  geom_point(aes(x=times, y=y.analytical, color="pink"))+
  xlab("Time")+
  ylab("X or Y")+
  theme_bw()

discrete.plot


# Continuous Time

# numerical solution

params <- c(a11, a12, a21, a22)

cont.time <- function(t, states, params){
  x=states[1]
  y=states[2]
  a11<-params[1]
  a12 <- params[2]
  a21<- params[3]
  a22 <-params[4]
  
  dxdt = a11*x+a12*y
  dydt = a21*x+a22*y
  
  return(list(c(dxdt,dydt)))
}

initial.values <- c(initial.x,initial.y)
times <- seq(0,10, by=0.01)

results <- lsoda(initial.values, times, cont.time, params)
colnames(results) <- c("time","x","y")

results.df <- data.frame(results)

results.plot <- ggplot(results.df)+
  geom_line(aes(x=time, y=x), color="blue")+
  geom_line(aes(x=time,y=y), color="red")+
  theme_bw()

results.plot

# analytical solution

dx_analytical <- 0.5*exp(1.5*times)-0.5*exp(0.5*times)
dy_analytical <- 0.5*exp(1.5*times)+0.5*exp(0.5*times)


results.df <- results.df |>
  mutate("dx.analytical"=dx_analytical, "dy.analytical"=dy_analytical)

results.plot <- ggplot(results.df)+
  geom_line(aes(x=time, y=x), color="blue")+
  geom_line(aes(x=time,y=y), color="red")+
  geom_line(aes(x=time, y=dx.analytical), color="green")+
  theme_bw()

results.plot
