#read the data
library(MASS)
lh = readMat("../data/LH.mat")

#mu is the mean of the initial fluorescence, and also its variance (up to a scaling factor)
m.t = abs(rnorm(n=1, mean=0, sd=0.01))
m.c = abs(rnorm(n=1, mean=0, sd=0.01))

#kappa is the scale of the precision of the initial fluorescence
k = rgamma(n=1, shape=0.1, rate=0.00001)

#alpha and beta are logistic regression parameters that model the reaction efficiency decreasing with cycle number
a.t = rnorm(n=1, mean=10, sd=7.5)
a.c = rnorm(n=1, mean=10, sd=7.5)
b.t = rnorm(n=1, mean=-1.5, sd=0.5)
b.c = rnorm(n=1, mean=-1.5, sd=0.5)
    
#tau.p is the precision of the logistic regression model for p
t.p = rgamma(n=1, shape=0.5, rate=0.0005)

#tau.y is the precision of the fluorescence observation
t.y = rgamma(n=1, shape=0.5, rate=0.0005)


#The total number of reactions run:
n.t = dim(lh[['TargetMatrix']])[2]
n.c = dim(lh[['StdMatrix']])[2]
n = n.t + n.c


#Draw the initial values of X0:
X0 = list()
X0[[1]] = list()
X0[[1]][['t']] = rnorm(n=n.t, mean=m.t, sd=sqrt(m.t/k))
X0[[1]][['c']] = rnorm(n=n.c, mean=m.c, sd=sqrt(m.c/k))


#Number of exponential-phase cycles per reaction:
first = 17
last = 20
m = last-first+1

#Start the chain with X=Y (noiseless observations)
X = Y = list()
X[[1]] = list()
Y[['t']] = X[[1]][['t']] = lh[['TargetMatrix']][first:last,]
Y[['c']] = X[[1]][['c']] = lh[['StdMatrix']][first:last,]

#Add the X0's to X
X[[1]][['t']] = rbind(X0[[1]][['t']], X[[1]][['t']])
X[[1]][['c']] = rbind(X0[[1]][['c']], X[[1]][['c']])

#Set our initial p
p = list()
p[[1]] = list()
#p.t = exp(a.t) / (exp(a.t) + X[[1]][['t']]**-b.t)
#p.c = exp(a.c) / (exp(a.c) + X[[1]][['c']]**-b.c)

p.t = matrix(0.99, nrow=4, ncol=15)
p.c = matrix(0.99, nrow=4, ncol=15)
p[[1]][['t']] = p.t
p[[1]][['c']] = p.c


#MCMC parameters
tau.y = vector()
tau.p = vector()
kappa = vector()
mu.t = vector()
mu.c = vector
alpha.t = vector()
beta.t = vector()
alpha.c = vector()
beta.c = vector()

mu.t = c(m.t)
mu.c = c(m.c)
alpha.t = c(a.t)
beta.t = c(b.t)
alpha.c = c(a.c)
beta.c = c(b.c)

#This function is used in updating tau.y
tauysum = function(Y, X, p) {
    tot = 0
    for (k in 1:dim(Y)[2]) {
        for (l in 1:dim(Y)[1]) {
            tot = tot + (Y[l,k] - X[k]*prod(1+p[1:l,k]))**2
        }
    }
    return(tot)
}
        
#This function is used in updating tau.p
taupsum = function(Y, X, p, a, b) {
    tot = 0
    for (j in 1:dim(Y)[2]) {
        for (i in 1:dim(Y)[1]) {
            tot = tot + (log(p[i,j]/(1-p[i,j])) - a - b*X[j]*prod(1+p[1:i,j]))**2
        }
    }
    return(tot)
}


#This function generates the X matrix for the regression parameters a, b
Xmat = function(X, p) {
    X2 = vector()
    for (i in 1:dim(p)[1]) {
        for (j in 1:dim(p)[2]) {
            X2 = c(X2, X[j]*prod(1+p[1:i,j]))
        }
    }
    return(as.matrix(data.frame(X1=rep(1,length(X2)), X2=X2)))
}

#This function generates the X matrix for the regression parameters a, b
getEta = function(p) {
    eta = vector()
    for (i in 1:dim(p)[1]) {
        for (j in 1:dim(p)[2]) {
            eta = c(eta, log(p[i,j]/(1-p[i,j])))
        }
    }
    return(eta)
}


#Get the exponent for the likelihood for the Metropolis step for finding the next p_ij
getPExp = function(X0, p, Y, a, b, i, j, t.y, t.p) {
    tot = 0
    m = dim(p)[1]
    for (k in i:m) {
        tot = tot + t.p/2 * (log(p[k,j]/(1-p[k,j])) - a - b*X0[j]*prod((1+p[1:k,j]))/(1+p[k,j]))**2
        tot = tot + t.y/2 * (Y[k,j] - X0[j]*prod(1+p[1:k,j]))**2
    }
    return(tot)
}


#Get the exponent for the likelihood for the Metropolis step for finding the next X0_j
getX0Exp = function(X0, p, Y, a, b, j, t.y, t.p) {
    tot = 0
    m = dim(p)[1]
    for (k in 1:m) {
        tot = tot + t.p/2 * (log(p[k,j]/(1-p[k,j])) - a - b*X0[j]*prod((1+p[1:k,j]))/(1+p[k,j]))**2
        tot = tot + t.y/2 * (Y[k,j] - X0[j]*prod(1+p[1:k,j]))**2
    }
    return(tot)
}


accept.u.t = vector()
accept.u.c = vector()
S = 100000

for (i in 2:S) {
    #Print the iteration number:
    print(i)

    #Draw new precision parameters:
    k.next = rgamma(n=1, shape=0.1+n/2, rate=0.00001+sum((X0[[i-1]][['c']]-mu.c[i-1])**2/(2*mu.c[i-1])) + sum((X0[[i-1]][['t']]-mu.t[i-1])**2/(2*mu.t[i-1])))
    t.y.next = rgamma(n=1, shape=0.5+m*n/2, rate=0.0005+0.5*(tauysum(Y[['c']], X0[[i-1]][['c']], p[[i-1]][['c']]) + tauysum(Y[['t']], X0[[i-1]][['t']], p[[i-1]][['t']])))
    t.p.next = rgamma(n=1, shape=0.5+m*n/2, rate=0.0005+0.5*(taupsum(Y[['c']], X0[[i-1]][['c']], p[[i-1]][['c']], a.c, b.c) + taupsum(Y[['t']], X0[[i-1]][['t']], p[[i-1]][['t']], a.t, b.t)))

    #propose a gaussian jump for mu.t:
    J = rnorm(n=1, mean=mu.t[i-1], sd=0.005)
    R = exp(-500*J**2 - k.next/(2*J)*sum((X0[[i-1]][['t']] - J)**2)) * mu.t[i-1] / (exp(-500*mu.t[i-1]**2 - k.next/(2*mu.t[i-1])*sum((X0[[i-1]][['t']] - mu.t[i-1])**2)) * J)
    U = runif(1)
    if (J<0) {m.t.next = mu.t[i-1]
    } else if (R>1) {m.t.next = J
    } else if (R>U) {m.t.next = J
    } else {m.t.next = mu.t[i-1]}

    #propose a gaussian jump for mu.c:
    J = rnorm(n=1, mean=mu.c[i-1], sd=0.005)
    R = exp(-500*J**2 - k.next/(2*J)*sum((X0[[i-1]][['c']] - J)**2)) * mu.c[i-1] / (exp(-500*mu.c[i-1]**2 - k.next/(2*mu.c[i-1])*sum((X0[[i-1]][['c']] - mu.c[i-1])**2)) * J)
    U = runif(1)
    if (J<0) { m.c.next = mu.c[[i-1]]
    } else if (R>1) { m.c.next = J 
    } else if (R>U) { m.c.next = J
    } else {m.c.next = mu.c[[i-1]]}

    #Logistic regression coefficients for the treatment group:
    mat = Xmat(X0[[i-1]][['t']], p[[i-1]][['t']])
    eta = getEta(p[[i-1]][['t']])
    M1 = as.matrix(lm(eta~mat-1)[['coefficients']])

    M2 = matrix(c(10, -1.5), 2, 1)
    P2 = matrix(c(7.5**-2, 0, 0, 0.5**-2),2,2)
    M = t(t(t.p.next * solve(t(mat) %*% mat) %*% M1 + P2 %*% M2) %*% solve(P2 + t.p.next * solve(t(mat) %*% mat)))
    
    ab.t.next = mvrnorm(n=1, mu=M, Sigma=solve(P2 + t.p.next * solve(t(mat) %*% mat)))
    a.t.next = ab.t.next[1]
    b.t.next = ab.t.next[2]


    #Logistic regression coefficients for the control group:
    mat = Xmat(X0[[i-1]][['c']], p[[i-1]][['c']])
    eta = getEta(p[[i-1]][['c']])
    M1 = as.matrix(lm(eta~mat-1)[['coefficients']])

    M2 = matrix(c(10, -1.5), 2, 1)
    P2 = matrix(c(7.5**-2, 0, 0, 0.5**-2),2,2)
    M = t(t(t.p.next * solve(t(mat) %*% mat) %*% M1 + P2 %*% M2) %*% solve(P2 + t.p.next * solve(t(mat) %*% mat)))
    
    ab.c.next = mvrnorm(n=1, mu=M, Sigma=solve(P2 + t.p.next * solve(t(mat) %*% mat)))
    a.c.next = ab.c.next[1]
    b.c.next = ab.c.next[2]



    #propose a gaussian jump for each p (control group):
    orig = p.c.next = p[[i-1]][['c']]
    for (j in 1:15) {
        for (k in 1:4) {
            J = rnorm(n=1, mean=p[[i-1]][['c']][k,j], sd=0.005)
            if(J>0 && J<1) {
                proposed = p.c.next
                proposed[k,j] = J
                R = exp(-getPExp(X0[[i-1]][['c']], proposed, Y[['c']], a.c.next, b.c.next, k, j, t.y.next, t.p.next)) * p[[i-1]][['c']][k,j] * (1-p[[i-1]][['c']][k,j]) / (exp(-getPExp(X0[[i-1]][['c']], p.c.next, Y[['c']], a.c.next, b.c.next, k, j, t.y.next, t.p.next)) * J * (1-J))
                U = runif(1)
                if (R>U) { p.c.next = proposed }
            }
        }
    }
    

    #propose a gaussian jump for each p (treatment group):
    orig = p.t.next = p[[i-1]][['t']]
    for (j in 1:15) {
        for (k in 1:4) {
            J = rnorm(n=1, mean=p[[i-1]][['t']][k,j], sd=0.005)
            if(J>0 && J<1) {
                proposed = p.t.next
                proposed[k,j] = J
                R = exp(-getPExp(X0[[i-1]][['t']], proposed, Y[['t']], a.t.next, b.t.next, k, j, t.y.next, t.p.next)) * p[[i-1]][['t']][k,j] * (1-p[[i-1]][['t']][k,j]) / (exp(-getPExp(X0[[i-1]][['t']], p.t.next, Y[['t']], a.t.next, b.t.next, k, j, t.y.next, t.p.next)) * J * (1-J))
                U = runif(1)
                if (R>U) { p.t.next = proposed }
            }
        }
    }





    #propose a gaussian jump for each X0 (control group):
    orig = X0.c.next = X0[[i-1]][['c']]
    for (j in 1:15) {
        J = rnorm(n=1, mean=orig[j], sd=0.005)
        proposed = X0.c.next
        proposed[j] = J
        R = exp(-getX0Exp(proposed, p.c.next, Y[['c']], a.c.next, b.c.next, j, t.y.next, t.p.next)) / exp(-getX0Exp(X0.c.next, p.c.next, Y[['c']], a.c.next, b.c.next, j, t.y.next, t.p.next))
        U = runif(1)
        if (R>U) { X0.c.next = proposed }
    }


    #propose a gaussian jump for each X0 (control group):
    orig = X0.t.next = X0[[i-1]][['t']]
    for (j in 1:15) {
        J = rnorm(n=1, mean=orig[j], sd=0.005)
        proposed = X0.t.next
        proposed[j] = J
        R = exp(-getX0Exp(proposed, p.t.next, Y[['t']], a.t.next, b.t.next, j, t.y.next, t.p.next)) / exp(-getX0Exp(X0.t.next, p.t.next, Y[['t']], a.t.next, b.t.next, j, t.y.next, t.p.next))
        U = runif(1)
        if (R>U) { X0.t.next = proposed }
    }


    kappa = c(kappa, k.next)
    tau.y = c(tau.y, t.y.next)
    tau.p = c(tau.p, t.p.next)
    mu.t = c(mu.t, m.t.next)
    mu.c = c(mu.c, m.c.next)
    alpha.t = c(alpha.t, a.t.next)
    beta.t = c(beta.t, b.t.next)
    alpha.c = c(alpha.c, a.c.next)
    beta.c = c(beta.c, b.c.next)
    p[[i]] = list()
    p[[i]][['c']] = p.c.next
    p[[i]][['t']] = p.t.next
    X0[[i]] = list()
    X0[[i]][['c']] = X0.c.next
    X0[[i]][['t']] = X0.t.next
    
    
}