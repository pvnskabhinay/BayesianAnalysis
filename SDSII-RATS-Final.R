# RATS - FINAL
library(rjags)
library(R2jags)
library(mcmcplots, quietly = T)
library(ggmcmc, quietly = T)
library(corrplot, quietly = T)

# I download the data from https://people.maths.bris.ac.uk/~mazjcr/BMB/2016/home.html
rats.data <- read.csv("D:\\Chrome Downloads\\SDS II\\Pochiraju\\Rats.csv", sep = ",",encoding = "utf-8")
head(rats.data)
# t1 = 8
# t2 = 15
# t3 = 22
# t4 = 29
# t5 = 36
t = c(8,15,22,29,36)
N = 30
# Visualize the data
matplot(t, t(rats.data), type = "b", pch = 16, lty = 1)
title(main = "Growth-rate of each Rat for 5-weeks")
grid() 

#rats.data.list = list(t = t, tbar = 22, N = 30, T = 5, Y = structure(.Data = rats.data,.dim=c(30,5)))
rats.data.raw <- list(x = c(8.0, 15.0, 22.0, 29.0, 36.0),
                  N = 30,
                  T = 5,
                  xbar=22,
                  Y = matrix(c(151, 199, 246, 283, 320,
                               145, 199, 249, 293, 354,
                               147, 214, 263, 312, 328,
                               155, 200, 237, 272, 297,
                               135, 188, 230, 280, 323,
                               159, 210, 252, 298, 331,
                               141, 189, 231, 275, 305,
                               159, 201, 248, 297, 338,
                               177, 236, 285, 350, 376,
                               134, 182, 220, 260, 296,
                               160, 208, 261, 313, 352,
                               143, 188, 220, 273, 314,
                               154, 200, 244, 289, 325,
                               171, 221, 270, 326, 358,
                               163, 216, 242, 281, 312,
                               160, 207, 248, 288, 324,
                               142, 187, 234, 280, 316,
                               156, 203, 243, 283, 317,
                               157, 212, 259, 307, 336,
                               152, 203, 246, 286, 321,
                               154, 205, 253, 298, 334,
                               139, 190, 225, 267, 302,
                               146, 191, 229, 272, 302,
                               157, 211, 250, 285, 323,
                               132, 185, 237, 286, 331,
                               160, 207, 257, 303, 345,
                               169, 216, 261, 295, 333,
                               157, 205, 248, 289, 316,
                               137, 180, 219, 258, 291,
                               153, 200, 244, 286, 324),
                             nrow=30, ncol=5, byrow=T))


Y <- rats.data.raw$Y
Y
x <-rats.data.raw$x
x
xbar <- rats.data.raw$xbar
xbar
N <-rats.data.raw$N
N
T <- rats.data.raw$T
T
# Defining the model

# 2.1 Frequentist Approach

c <- rep(NA,N) #intercept
m <- rep(NA,N) #slope of the model
for(i in 1:N)
{
  lm.Frequentist <- lm(rats.data.raw$Y[i,] ~ rats.data.raw$x)
  c[i] <- lm.Frequentist$coefficients[[1]]
  m[i] <- lm.Frequentist$coefficients[[2]]
}
c
m

#mean of the intercept
mean.alpha <- mean(c)
mean.alpha
#106.5676

#mean of the slope
mean.beta <- mean(m)
mean.beta
#6.185714

plot(rats.data.raw$x,colMeans(rats.data.raw$Y), lwd=4, xlab = "age(day)", ylab = "weight",
     col="red", ylim=c(135,355))
points(rep(rats.data.raw$x[1],N), rats.data.raw$Y[,1])
points(rep(rats.data.raw$x[2],N), rats.data.raw$Y[,2])
points(rep(rats.data.raw$x[3],N), rats.data.raw$Y[,3])
points(rep(rats.data.raw$x[4],N), rats.data.raw$Y[,4])
points(rep(rats.data.raw$x[5],N), rats.data.raw$Y[,5])
abline(mean.alpha, mean.beta, col="purple", lwd=2)
grid()


# 2.2 Normal Hierarchial model with linear expected value

model1<- function()
{
  for (i in 1:N)
  {
    for (j in 1:T)
    {
      Y[i,j] ~ dnorm(mu[i,j], tau.c)
      mu[i, j] <- alpha[i] + beta[i] * (x[j])
    }
    alpha[i] ~ dnorm(alpha.c, tau.alpha)
    beta[i] ~ dnorm(beta.c, tau.beta)
  }
  #Priors
  alpha.c ~ dnorm(0, 1.0E-6)
  beta.c ~ dnorm(0, 1.0E-6)
  tau.c ~ dgamma(1.0E-3, 1.0E-3)
  tau.alpha ~ dgamma(1.0E-3, 1.0E-3)
  tau.beta ~ dgamma(1.0E-3, 1.0E-3)
  #Transformation
  sigma.c <- 1.0/sqrt(tau.c)
  xbar <- mean(x[])
  alpha0 <- alpha.c - beta.c*xbar
}

#read rats data for Jags
rats.data.jags.list <- list("Y", "x", "T", "N")

# JAGS parameters
rats.params <- c("tau.c", "alpha.c", "beta.c", "tau.alpha", "tau.beta")

## Define the starting values for JAGS
rats.inits.1 <- function(){
  list("alpha"=c(250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250), 
       "beta"=c(6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6), 
       "alpha.c"= 150, "beta.c"= 10, "tau.c"= 1, "tau.alpha"= 1, "tau.beta"= 1)}

## Fit the model in JAGS, having previously copied the BUGS model in my working directory as "rats.model.jags"
ratsfit.model1 <- jags(data=rats.data.jags.list, inits=rats.inits.1, rats.params, n.chains=2, n.iter=10000,
                       n.burnin=1000, n.thin = 1, model.file=model1, DIC=TRUE)

ratsfit.model1
summary(ratsfit.model1)

ratsfit.model1$BUGSoutput$DIC #lower is good

#now write the mean of theintercept for this model
ratsfit.model1$BUGSoutput$summary[,"mean"]["alpha.c"]

#intercept
ratsfit.model1$BUGSoutput$summary[,"mean"]["beta.c"]

#plot the model1 just like frequentist model
plot(rats.data.raw$x,colMeans(rats.data.raw$Y), lwd=4, xlab = "age(days)", ylab = "weight",
     col="red", ylim=c(135,355))
points(rep(rats.data.raw$x[1],N), rats.data.raw$Y[,1])
points(rep(rats.data.raw$x[2],N), rats.data.raw$Y[,2])
points(rep(rats.data.raw$x[3],N), rats.data.raw$Y[,3])
points(rep(rats.data.raw$x[4],N), rats.data.raw$Y[,4])
points(rep(rats.data.raw$x[5],N), rats.data.raw$Y[,5])
abline(ratsfit.model1$BUGSoutput$summary[,"mean"]["alpha.c"],
       ratsfit.model1$BUGSoutput$summary[,"mean"]["beta.c"], col="orange", lwd=2)
grid()

#plot the model to see all the values
plot(ratsfit.model1)

traceplot(ratsfit.model1)

# 3. A normal Hierarchial model with different priors tau.alpha and tau.beta
model2 <- function()
{
  for (i in 1:N)
  {
    for (j in 1:T)
    {
      Y[i,j] ~ dnorm(mu[i,j], tau.c)
      mu[i, j] <- alpha[i] + beta[i] * (x[j])
    }
    alpha[i] ~ dnorm(alpha.c, tau.alpha)
    beta[i] ~ dnorm(beta.c, tau.beta)
  }
  alpha.c ~ dnorm(0, 1.0E-6)
  beta.c ~ dnorm(0, 1.0E-6)
  tau.c ~ dgamma(1.0E-3, 1.0E-3)
  sigma.alpha ~ dunif(0,100)
  sigma.beta ~ dunif(0,100)
  tau.alpha <- 1/(sigma.alpha*sigma.alpha)
  tau.beta <- 1/(sigma.beta*sigma.beta)
  sigma.c <- 1.0/sqrt(tau.c)
  xbar <- mean(x[])
  alpha0 <- alpha.c - beta.c*xbar
}
## Read in the rats data for JAGS
rats.data.list <- list("Y", "x", "T", "N")
## Name the JAGS parameters
rats.params <- c("tau.c", "alpha.c", "beta.c", "tau.alpha", "tau.beta")
#### Define the starting values for JAGS
rats.inits.2 <- function(){
  list(alpha = c(250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250,
                 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250),
       beta = c(6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6),
       alpha.c = 150, beta.c = 10,
       tau.c = 1, sigma.alpha = 1, sigma.beta = 1)
}
ratsfit.model2 <- jags(data=rats.data.list, inits=rats.inits.2, rats.params, n.chains=2, n.iter=10000,
             n.burnin=1000, n.thin = 1, model.file=model2, DIC=TRUE)

ratsfit.model2

ratsfit.model2$BUGSoutput$DIC #lower is good

#now write the mean of theintercept for this model
ratsfit.model2$BUGSoutput$summary[,"mean"]["alpha.c"]

#intercept
ratsfit.model2$BUGSoutput$summary[,"mean"]["beta.c"]

#plot the model1 just like frequentist model

plot(rats.data.raw$x,colMeans(rats.data.raw$Y), lwd=4, xlab = "age(days)", ylab = "weight",
     col="red", ylim=c(135,355))
points(rep(rats.data.raw$x[1],N), rats.data.raw$Y[,1])
points(rep(rats.data.raw$x[2],N), rats.data.raw$Y[,2])
points(rep(rats.data.raw$x[3],N), rats.data.raw$Y[,3])
points(rep(rats.data.raw$x[4],N), rats.data.raw$Y[,4])
points(rep(rats.data.raw$x[5],N), rats.data.raw$Y[,5])
abline(ratsfit.model2$BUGSoutput$summary[,"mean"]["alpha.c"],
       ratsfit.model2$BUGSoutput$summary[,"mean"]["beta.c"], col="black", lwd=2)
grid()
#plot the model to see all the values
plot(ratsfit.model2)


xyplot(ratsfit.model2)
traceplot(ratsfit.model2)


densityplot(ratsfit.model2)


#now compare the model with DIC
# Compare the 2 models using DIC
DIC.Models <- cbind(model1.DIC=ratsfit.model1$BUGSoutput$DIC,
                    model2.DIC=ratsfit.model2$BUGSoutput$DIC)

DIC.Models

#compare dic with plots

DIC = c(ratsfit.model1$BUGSoutput$DIC, ratsfit.model2$BUGSoutput$DIC)

barplot(c(ratsfit.model1$BUGSoutput$DIC, ratsfit.model2$BUGSoutput$DIC), col=c("blue", "red"), main="DIC comparison of Model 1 & 2",
        names.arg = c("model1", "model2"), ylim=c(0,1250))



plot(rats.data.raw$x,colMeans(rats.data.raw$Y), lwd=4, xlab = "age(day)", ylab = "weight",
     col="red", ylim=c(135,355))
points(rep(rats.data.raw$x[1],N), rats.data.raw$Y[,1])
points(rep(rats.data.raw$x[2],N), rats.data.raw$Y[,2])
points(rep(rats.data.raw$x[3],N), rats.data.raw$Y[,3])
points(rep(rats.data.raw$x[4],N), rats.data.raw$Y[,4])
points(rep(rats.data.raw$x[5],N), rats.data.raw$Y[,5])
abline(mean.alpha, mean.beta, col="purple", lwd=2)
abline(ratsfit.model1$BUGSoutput$summary[,"mean"]["alpha.c"],
       ratsfit.model1$BUGSoutput$summary[,"mean"]["beta.c"], col="orange", lwd=2)
abline(ratsfit.model2$BUGSoutput$summary[,"mean"]["alpha.c"],
       ratsfit.model2$BUGSoutput$summary[,"mean"]["beta.c"], col="red", lwd=2)
grid()


##see the convergence test
#check the convergence test for every parameter

# Convergence tests
#tau.c
par(mfrow=c(1,1))

tau.c.vals <- ratsfit.model1$BUGSoutput$sims.array[,1,"tau.c"]
plot(tau.c.vals, xlab = "iterations", main="trace plot for tau.c",type="l")

par(mfrow=c(1,2))
hist(tau.c.vals, main= "tau.c histogram", xlab = "tau.c")
abline(v=mean(tau.c.vals), col="red", lwd=2)
plot(cumsum(tau.c.vals)/(1:length(tau.c.vals)), type="l", ylab="",
     main="Convergence", xlab="simulations")
abline(h=mean(tau.c.vals), col="red", lwd=2)


#alpha.c
par(mfrow=c(1,1))

alpha.c.vals <- ratsfit.model1$BUGSoutput$sims.array[,1,"alpha.c"]
plot(alpha.c.vals, xlab = "iterations", main="alpha.c trace plot",type="l")
par(mfrow=c(1,2))
hist(alpha.c.vals, main= "alpha.c histogram", xlab = "alpha.c")
abline(v=mean(alpha.c.vals), col="red", lwd=2)
plot(cumsum(alpha.c.vals)/(1:length(alpha.c.vals)), type="l", ylab="",
     main="Convergence", xlab="simulations")
abline(h=mean(alpha.c.vals), col="red", lwd=2)

#beta.c
beta.c.vals <- ratsfit.model1$BUGSoutput$sims.array[,1,"beta.c"]

plot(beta.c.vals, xlab = "iterations", main="beta.c trace plot",type="l")
par(mfrow=c(1,2))

hist(beta.c.vals, main= "beta.c histogram", xlab = "beta.c")
abline(v=mean(beta.c.vals), col="red", lwd=2)

plot(cumsum(beta.c.vals)/(1:length(beta.c.vals)), type="l", ylab="",
     main="Convergence", xlab="simulations")
abline(h=mean(beta.c.vals), col="red", lwd=2)

par(mfrow=c(1,1))


#tau.alpha
tau.alpha.vals <- ratsfit.model1$BUGSoutput$sims.array[,1,"tau.alpha"]
plot(tau.alpha.vals, xlab = "iterations", main="tau.alpha trace plot",type="l")

par(mfrow=c(1,2))
hist(tau.alpha.vals, main= "tau.alpha histogram", xlab = "tau.alpha")
abline(v=mean(tau.alpha.vals), col="red", lwd=2)

plot(cumsum(tau.alpha.vals)/(1:length(tau.alpha.vals)), type="l", ylab="",
     main="Convergence", xlab="simulations")
abline(h=mean(tau.alpha.vals), col="red", lwd=2)

#tau.beta
par(mfrow=c(1,1))
# tau.beta
tau.beta.vals <- ratsfit.model1$BUGSoutput$sims.array[,1,"tau.beta"]
plot(tau.beta.vals, xlab = "iterations", main="tau.beta trace plot",type="l")

par(mfrow=c(1,2))
hist(tau.beta.vals, main= "tau.beta histogram", xlab = "tau.beta")
abline(v=mean(tau.beta.vals), col="red", lwd=2)

plot(cumsum(tau.beta.vals)/(1:length(tau.beta.vals)), type="l", ylab="",
     main="Convergence", xlab="simulations")
abline(h=mean(tau.beta.vals), col="red", lwd=2)

# Prediction
prediction <- function(x){
  alpha <- rep(NA, 9000)
  beta <- rep(NA, 9000)
  y.pred <- rep(NA, 9000)
  for(i in 1:9000){
    alpha[i] = rnorm(1,alpha.c.vals[i], tau.alpha.vals[i])
    beta[i] = rnorm(1,beta.c.vals[i], tau.beta.vals[i])
    y.pred[i] = alpha[i]+beta[i]*x
  }
  return(y.pred)
}

mean(prediction(0))
mean(prediction(8))
mean(prediction(15))
mean(prediction(22))
mean(prediction(29))
mean(prediction(36))





