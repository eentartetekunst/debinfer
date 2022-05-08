# install.packages("deSolve"); library("deSolve") 
# install.packages("DEoptim"); library("DEoptim")

myData <- read.csv(file = '/Users/sophielohshvits/Desktop/myData.csv')
myData

# UIV model of viral dynamics in the cell
myModel <- function(t,state,parameters) {
  with(as.list(c(state,parameters)),{
    dU = - myBeta*U*V
    dI = myBeta*U*V - myDelta*I 
    dV = myP*I - myC*V
    list(c(dU, dI, dV))
    })
}


myStates <- c(U = 10^6, I = 0, V = 10) # initial conditions
myParams <- c("myBeta","myDelta","myP","myC") # list of parameters 
myV <- "V" # the component in the model observed in data


modelTime <- seq(from = 0, to = 10, by = 0.01) # observing fata every 15 mins a day


# cost matrix.  The root mean square errors of the model output and experimental data
myCostFn <- function(x) {
  parms <- x[1:length(myParams)] 
  names(parms) <- myParams
  yhat <- ode(myStates, modelTime, myModel, 10^parms)
  yMatch <- yhat[yhat[,1] %in% myData$time, ] 
  nm <- rle(myData$time)$lengths
  x <- myData[, myV]- rep(log10(yMatch[, myV]), times = nm)
  rmse <- sqrt(mean(x^2))
  return(rmse) 
} 

# parameter boundaries and optimising conditions: 
lower = log10(c(1e-7, 1e-2, 1e+0, 1e-1)) # 
upper = log10(c(1e-3, 1e+2, 1e+2, 1e+2))

# half life of infected cell -- 
myOptions <- DEoptim.control(itermax = 10000, steptol = 50, reltol = 1e-8) # see Notes 12

fit <- do.call("DEoptim", list(myCostFn, lower, upper, myOptions))


(bestPar <- fit$optim$bestmem) # parameter estimates in log 10 scale 
names(bestPar) <- myParams
out <- ode(myStates, modelTime, myModel, 10^bestPar)
plot(out[, "time"], log10(out[, "V"]), type = "l") # in log 10
points(myData$time, myData$V) # superimposing observed data points


myAIC <- function(fit, np=NULL, rms=NULL, n=NULL) {
  if (is.null(n)) stop("How many observations were used? n=#") 
  if (is.null(np)) np <- length(fit$optim$bestmem)
  if (is.null(rms)) rms <- fit$optim$bestval
  return(2*np + n*log(rms))
}

myAIC(fit, n = 30)


myProfile <- function(lower, upper, bestPar) { 
  pro.ll <- NULL
  for (v in 1:length(bestPar)) {
    # Creating parameter sequence
    tmpl <- seq(lower[v], bestPar[[v]], length.out = 100)
    tmpl <- tmpl[order(tmpl, decreasing = TRUE)[cumsum(1:13)]] 
    tmpr <- seq(bestPar[[v]], upper[v], length.out = 100)
    tmpr <- tmpr[cumsum(1:13)]
    pars <- sort(unique(c(lower[v], tmpl, bestPar[[v]], tmpr, upper[v]))) 
    ppl <- NULL
    # Run optimization for each and record the parameters and RMSE
    for (p in pars) {
      DEargs <- list(myCostFn, replace(lower, v, p), replace(upper, v, p), myOptions)
      fit <- do.call("DEoptim", DEargs)
      ppl <- c(ppl, fit$optim$bestval) 
    }
    pro.ll[[v]] <- cbind(pars, ppl)
  } 
  return(pro.ll)
}


outProfiles <- myProfile(lower, upper, bestPar)

par(mfrow = c(2, 2))
sapply(1:4, function(x) plot(outProfiles[[x]], xlab = myParams[x], ylab = 'RMSE'))

myBoot <- function(numboot = 1000, numpar = 4) {
  # numpar: number of parameters in the model
  # numboot: number of bootstrap samples
  results <- matrix(NA, numboot, numpar) 
  original <- myData
  sampling <- function(x) sample(original$V[original$time==x], length(original$V[original$time==x]), replace = 1)
  for (i in 1:numboot) {
    message("Bootstrapping sample ", i)
    tmp <- sapply(unique(original$time), sampling) 
    myData <- cbind(original$time, as.vector(tmp)) 
    DEarguments <- list(myCostFn, lower, upper, myOptions) 
    fit <- do.call("DEoptim", DEarguments)
    results[i, ] <- fit$optim$bestmem
  }
  results <- as.data.frame(results) 
  colnames(results) <- myParams
  return(results) 
}

bootResults <- myBoot()

par(mfrow = c(2,2))
sapply(1:4, function(x) hist(bootResults[, x]) )

apply(bootResults, 2, quantile, probs = c(.025,.975))

