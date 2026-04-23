library(deSolve)
# Continuous-time linear constant coefficient ODEs 
a11 = -7
a12 = 15
a21 = -6
a22 = 11

params = c(a11, a12, a21, a22)
A<-matrix(params,nrow=2,ncol=2,byrow=TRUE)
A
eigen(A)


Conttime <- function (t, states, params) {
  X = states[1]
  Y = states[2]
  a11 = params[1]
  a12 = params[2]
  a21 = params[3]
  a22 = params[4]
  
  dXdt = a11*X + a12*Y
  dYdt = a21*X + a22*Y 
  return(list(c(dXdt,dYdt)))
}


X0 = 0 # Initial X
Y0 = 0.1 # Initial Y
initial_values = c(X=X0,Y=Y0)
params = c(a11, a12, a21, a22)

times = seq(0, 6, by=0.01)

results = lsoda(initial_values, times, Conttime, params)
colnames(results) = c("time", "X", "Y")

plot(results[, "time"], results[, "X"], type="l", col="blue", xlab="time",ylab="X or Y")
lines(results[, "time"], results[, "Y"], type="l", col="red")

############################
# Example from last Thursday
# dx/dt = 3x - y
# dy/dt = -2x + 2 y
# x(0) = 90
# y(0) = 150
############################

a11 = 3
a12 = -1
a21 = -2
a22 = 2

params = c(a11, a12, a21, a22)
A<-matrix(params,nrow=2,ncol=2,byrow=TRUE)
A
eigen(A)

#  Find solution numerically

X0 = 90 # Initial X
Y0 = 150 # Initial Y

initial_values = c(X=X0,Y=Y0)
params = c(a11, a12, a21, a22)

times = seq(0, 1, by=0.01)

results = lsoda(initial_values, times, Conttime, params)
colnames(results) = c("time", "X", "Y")

plot(results[, "time"], results[, "X"], type="l", col="blue", xlab="time",ylab="x or y",ylim=c(0,max(results[, "X"])))
lines(results[, "time"], results[, "Y"], type="l", col="red")

# We found the solution Analytically to be: 
# x(t) = 10*exp(4*t) + 80*exp(t)
# y(t) = -10*exp(4*t) + 160*exp(t)

# To plot the analytical solution, we don't need an ode solver. 
# We can just calculate it.
# We generated the vector of times above, so let's use that.

x = 10*exp(4*times) + 80*exp(times)
y = -10*exp(4*times) + 160*exp(times)

points(times, x, pch = 5, col = "blue")
points(times, y, pch = 1, col = "red")

legend("topleft", legend=c("X numerical", "X analytical", "Y numerical", "Y analytical"), col=c("blue", "blue", "red", "red"), , lty = c(1, NA, 1, NA), pch = c(NA, 5, NA, 1))


#################  Discrete time ##################
# Let's find numerical and analytical solutions to:
# x(t+1) = 3 x(t) - y(t)
# y(t+1) = -2 x(t) + 2 y(t)
# x(0) = 90
# y(0) = 150
####################################################

a11 = 3
a12 = -1
a21 = -2
a22 = 2

params = c(a11, a12, a21, a22)
M<-matrix(params,nrow=2,ncol=2,byrow=TRUE)
M
eigen(M)

# Discrete-time, numerical solution

# Need a for loop to get the numerical solution

initial_x = 90
initial_y = 150

max_time = 5

x = rep(0,(max_time+1))
y = rep(0,(max_time+1))
x[1] = initial_x
y[1] = initial_y
for (t in 1:max_time) {
  x[t+1] = a11*x[t] + a12*y[t]
  y[t+1] = a21*x[t] + a22*y[t]
}

times<-seq(0,max_time)
plot(times,x,type="l",xlab="time",ylab="X or Y",col="blue")
lines(times,y,col="red")


#Discrete-time, analytical solution

x_analytical = 10*4^times + 80*1^times
y_analytical = -10*4^times + 160*1^times

points(times, x_analytical, pch = 5, col = "blue")
points(times, y_analytical, pch = 1, col = "red")

legend("topleft", legend=c("X numerical", "X analytical", "Y numerical", "Y analytical"), col=c("blue", "blue", "red", "red"), lty = c(1, NA, 1, NA), pch = c(NA, 5, NA, 1))






