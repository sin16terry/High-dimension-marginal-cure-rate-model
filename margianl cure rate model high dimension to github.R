library(optimization)
library(MASS)
library(glmnet)
library(doParallel)
library(quad)
library(pracma) ##for incomplete gamma function gammainc
library(maxLik)

####################################################################################
################################set true parameters#################################
rm(list = ls(all = TRUE))
MC = 500          ## replicating times of simulation
mse = NA          ## mean squre error
n = 200           ## sample size
p = 300           ## number of covariates including the intercept term
salpha = 1.1      ## nusiance parameter alpha in the Weibull distribution
exp.theta = 6     ## nusiance parameter theta in the two-parameter exponential distribution
cen.rate = 0.0002     ##parameter in the exponential distr. for generatingcensoring time
False.Post.Ra.1 = NA  ##false positive rate then threshhold value =0.1
False.Nega.Ra.1 = NA  ##false negative rate then threshhold value =0.1
False.Post.Ra.05 = NA ##false positive rate then threshhold value =0.05
False.Nega.Ra.05 = NA ##false positive rate then threshhold value =0.05
ave.uncensor = NA     ##average uncensor rate for each sample
para.estimation = matrix(rep(0, MC * (p + 2)), ncol = MC, nrow = p + 2)  #define the parameter vector
t1 = Sys.time()  ##measure the runtime
m = 1            ## replication start with 1
breakk = FALSE   ## if certain error appear, break the loop
true.para = c(-1, 1, -1.5, 2, 1.75, -1.25, rep(0, p - 6))  #define true value of paramteres
############define variance-covariance matrix for covariates matrix ###########
sig = matrix(0, p, p)
for (i in 1:p)
{
  for (j in 1:p)
  {
    sig[i, j] = 0.2 ^ abs(i - j)
  }
}
###############start the replication of simulation with MC times ##############
while (m <= MC) {
  #############################################################################
  ############################PART1:generate data##############################
  trsize = n                              ###sample size
  trainx = mvrnorm(trsize, rep(0, p), sig)###covariate matrix
  cccc = trainx[, -1]              
  intecept = c(rep(1, trsize))    #####constant intercept case
  trainx = as.matrix(cbind(intecept, cccc)) ##### add exposure treatment to the covariate
  for (i in 1:p) {
    assign(paste("x", i, sep = ""), (trainx[, i]))
  }
  gm0 = rexp(trsize, exp.theta)
  gm = gm0 + 1  ####gm is two parameter exponential distribution
  uncure.rate = gm ^ (-1 / salpha)  ###uncure rate
  #true marginal mean hazard rate marginh
  marmeanh = exp(
    true.para[1] * x1 + true.para[2] * x2 + true.para[3] * x3 + true.para[4] *
      x4 + true.para[5] * x5 + true.para[6] * x6
  )
  condmean.lameu = (marmeanh / (salpha * uncure.rate * gamma(2 - 1 / salpha))) ^ salpha
  ss = runif(trsize) ## use inverse of CDF to generate the random values for the noncured
  #summary(marmeanh)
  uncure.t = (-log(1 - ss) / condmean.lameu) ^ (1 / salpha)  ##surviva time for noncured individuals
  #which.max(uncure.t)#uncure.rate[50]#summary(uncure.t)#hist(uncure.t)
  g = NA   ##cure or noncur indicator  g=1 noncure, g=0 cure
  for (ii in 1:trsize) {
    g[ii] = rbinom(1, 1, uncure.rate[ii])
  }
  mean(g)
  expcensor.t = rexp(trsize, cen.rate)
  censor = NA    ## censoring  or follow-up time for each individual
  for (i in 1:trsize) {
    censor[i] = min(300, expcensor.t[i])
  }
  cure.t = rep(10000000, trsize)
  truet = NA    ## true survival time for each individual
  truet[g == 1] = uncure.t[g == 1]
  truet[g == 0] = cure.t[g == 0]
  d = NA                          #### censoring indicator, d=1 then noncensored d=0 censored
  for (i in 1:trsize) {
    if (truet[i] <= censor[i]) {
      d[i] = 1
    } else {
      d[i] = 0
    }
  }
  obt = NA      
  
  for (i in 1:trsize) {
    if (d[i] == 1) {
      obt[i] = truet[i]
    } else {
      obt[i] = censor[i]
    }
  }
  
  #summary(d)
  t = obt   ###observed time for each individual
  #hist(obt)
  ave.uncensor[m] = mean(d)
  ##########################  end data generating ############################
  ############################################################################
  
  
  
  #############################################################################
  ##################### PART2(maxmize the loglikelihood)#######################
  CONTINUE = TRUE  #### if not \beta converge, then false
  ITER = 0
  update = c(
    salpha,
    exp.theta,
    true.para[1],
    true.para[2],
    true.para[3],
    true.para[4],
    true.para[5],
    true.para[6],
    rep(0, p - 6)
  )
  Vgamma = Vectorize(gammainc)  ##### Vectorize the incomplete gamma function
  ####OUTERLOOP
  while (CONTINUE) {
    ITER = ITER + 1
    beta1 = update
    update.nuisance.ll = function (para) {
      #para=opti.nuisance$par
      #para=c(salpha,exp.theta)
      alpha = para[1]
      theta = para[2]
      sumeta = trainx %*% beta1[-(1:2)]
      ##summary(sumeta)
      mm = exp(sumeta)
      lameu = (mm / (alpha * gamma(2 - (1 / alpha)))) ^ alpha
      b = (t ^ alpha) * lameu
      ####find log likelihood when d==1
      ll = NA
      c = NA
      c[b <= 600] = b[b <= 600]
      c[b > 600] = 600
      ##################using  incomplete gamma since 2-1/alpha>0###############
      l2 = Vgamma(c + theta, 2 - 1 / alpha)[2, ]
      ll1 = log(alpha) + log(c) - log(t) + log(theta) + theta + (1 / alpha -
                                                                   2) * log(theta + b) + log(l2)
      ll[d == 1] = ll1[d == 1]
      sum(ll[d == 1])
      ####find log likelihood when d==0
      ##################using integrate function to calculation incomplete gamma#######
      ##################since 1-1/alpha<0 when alpha<1     ############################
      incom.gamma = function(x) {
        x ^ (-1 / alpha) * exp(-x)
      }
      inte.1pi = integrate(incom.gamma, 1, 700)$value - integrate(incom.gamma, 1, theta)$value
      e1 = 1 - theta ^ (1 / alpha) * exp(theta) * inte.1pi
      e2 = NA
      inte.pis = NA
      for (i in 1:trsize) {
        inte.pis[i] = integrate(incom.gamma, 1, 700)$value - integrate(incom.gamma, 1, c[i] +
                                                                         theta)$value
        e2[i] = (b[i] + theta) ^ (1 / alpha - 1) * theta * exp(theta) * inte.pis[i]
      }
      ll[d == 0] = log((e1 + e2)[d == 0])
      sum(ll[d == 0])
      sum(ll)
    }
    
    ####update.nuisance.ll(c(slambda,salpha,exp.theta))  ##TRUE NUISANCE PARAMETER
    ini = beta1[1:2]
    A <- matrix(c(1, 0, 0, 1), 2, 2)
    B <- c(-0.51, 0)
    opti.nuisance <-
      maxLik(
        update.nuisance.ll,
        start = ini,
        constraints = list(ineqA = A, ineqB = B)
      )
    update.nuisance.ll(opti.nuisance$estimate) - update.nuisance.ll(c(salpha, exp.theta))
    nuisance.update = opti.nuisance$estimate
    ###############End the estimation of nuisance parameter alpha & theta################
    #####################################################################################
    
    
    #####################################################################################
    ###################PART3(start the quadratic approximation)##########################
    #####################################################################################
    continue = TRUE   #
    iter = 0
    highdi.update = c(nuisance.update[1], nuisance.update[2], update[3:length(update)])
    ##### inner update the highdimension parameters beta
    while (continue) {
      iter = iter + 1
      beta = highdi.update
      alpha = beta[1]
      theta = beta[2]
      #create likelihood function input is eta
      llder = function(eta) {
        mm = exp(eta)
        lameu = (mm / (alpha * gamma(2 - 1 / alpha))) ^ alpha
        b = (t ^ alpha) * lameu
        ####find log likelihood when d==1
        ll = NA
        c = NA
        c[b <= 600] = b[b <= 600]
        c[b > 600] = 600
        l2 = Vgamma(c + theta, 2 - 1 / alpha)[2, ]
        ll1 = log(alpha) + log(c) - log(t) + log(theta) + theta + (1 / alpha -
                                                                     2) * log(theta + b) + log(l2)
        ll[d == 1] = ll1[d == 1]
        ####find log likelihood when d==0
        incom.gamma = function(x) {
          x ^ (-1 / alpha) * exp(-x)
        }
        inte.1pi = integrate(incom.gamma, 1, 700)$value - integrate(incom.gamma, 1, theta)$value
        e1 = 1 - theta ^ (1 / alpha) * exp(theta) * inte.1pi
        e2 = NA
        inte.pis = NA
        for (i in 1:trsize) {
          inte.pis[i] = integrate(incom.gamma, 1, 700)$value - integrate(incom.gamma, 1, c[i] +
                                                                           theta)$value
          e2[i] = (b[i] + theta) ^ (1 / alpha - 1) * theta * exp(theta) * inte.pis[i]
        }
        ll[d == 0] = log((e1 + e2)[d == 0])
        ll
      }
      #### calculate the first deravative wrt eta on xb
      eta = trainx %*% beta[-(1:2)]
      sum(llder(eta))
      h = 0.001
      mu = (llder(eta + h) - llder(eta - h)) / (2 * h)
      secder2 = -((llder(eta + h) + llder(eta - h) - 2 * llder(eta)) / (h ^
                                                                          2))
      breakk = any(is.nan(secder2))
      if (breakk) {
        break
      }
      #### delete points where second derivative <0
      summary(secder2)
      secder = NA
      secder[secder2 <= 0] = 0
      secder[secder2 > 0] = secder2[secder2 > 0]
      A = matrix(0, trsize, trsize)
      for (i in 1:trsize) {
        A[i, i] = secder[i]
      }
      
      #A[is.nan(A)]=0
      Q = A ^ 0.5
      z = eta + ginv(A) %*% mu
      ztuta = Q %*% z
      Xtuta = Q %*% trainx
      
      lasso.cv <-
        cv.glmnet(
          Xtuta,
          ztuta,
          intercept = FALSE,
          alpha = 1,
          nfolds = 5,
          type.measure = "deviance"
        )
      minalpha <- lasso.cv$lambda.min  # lambda in the notes
      lasso.fit <-
        glmnet(
          Xtuta,
          ztuta,
          family = "gaussian",
          intercept = FALSE,
          lambda = minalpha,
          alpha = 1
        )
      coef0 = coef(lasso.fit)[-1]
      #coef0 <- as.vector(predict(lasso.fit, s=minalpha,intercept=FALSE, type="coefficients"))[-1]
      highdi.update = c(beta[1], beta[2], coef0)
      res = highdi.update - beta
      continue = (crossprod(res, res) > 0.00001) & (iter <= 100)
      #print(iter)
    }
    if (breakk) {
      break
    }
    update = highdi.update
    res2 = update - beta1
    CONTINUE = (crossprod(res2, res2) > 0.00001) & (ITER <= 17)
    print(paste("current MCsample is", m, ITER))
  }
  if (breakk) {
    next
  }
  mse[m] = crossprod(update[-c(1, 2)] - true.para, update[-c(1, 2)] - true.para)
  para.estimation[, m] = update
  #####0.1 as cutting for calculating false positive
  final.para.1 = NA
  final.para.1[abs(update) < 0.1] = 0
  final.para.1[abs(update) > 0.1] = update[abs(update) > 0.1]
  False.Post.Ra.1[m] = (p - 6 - sum(final.para.1[-c(1:8)] == 0)) / (p -
                                                                      6)
  False.Nega.Ra.1[m] = sum(final.para.1[c(3:8)] == 0) / 6
  #####0.05 as cutting for calculating false positive
  final.para.05 = NA
  final.para.05[abs(update) < 0.05] = 0
  final.para.05[abs(update) > 0.05] = update[abs(update) > 0.05]
  False.Post.Ra.05[m] = (p - 6 - sum(final.para.05[-c(1:8)] == 0)) / (p -
                                                                        6)
  False.Nega.Ra.05[m] = sum(final.para.05[c(3:8)] == 0) / 6
  print(Sys.time() - t1)
  m = m + 1
}
std.est = apply(para.estimation, 1, std)
matrix.count = para.estimation
matrix.count[abs(para.estimation) < 0.05] = 0
matrix.count[abs(para.estimation) > 0.05] = 1
count.posi = apply(matrix.count, 1, sum)
mean(ave.uncensor)
mean(mse)
mean(std.est[c(3:8)])
mean(False.Post.Ra.1)
mean(False.Nega.Ra.1)
mean(False.Post.Ra.05)
mean(False.Nega.Ra.05)
print(para.estimation)
#write.csv(para.estimation,'healthstudy2.csv')
#fix(para.estimation)
print(count.posi)

###############################END SIMULATION#############################
##########################################################################
