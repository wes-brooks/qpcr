#Define the detectability threshold:
threshold = 10000

#Define the actual parameters of the process:
#Define the intensity of the processes that generate particle counts (we'd like to recover these):
lambda = list()
lambda[[1]] = 10
lambda[[2]] = 15

#Actual parameters for distribution of replication probabilities:
p.a = list()
p.a[[1]] = 96
p.a[[2]] = 93

#Reaction-level random effect parameter for efficiency and "sample size" for beta distribution:
p.sigma = 2
beta.ss = 100

#Number of reactions to be run on each gene/treatment combo:
R = 3

#Number of replication cycles per reaction:
N = 25

#Additive noise:
x.sigma = 5000

#Generate random effects for efficiency:
alpha = list()
beta = list()

alpha[[1]] = p.a[[1]] + rnorm(R, 0, p.sigma)
alpha[[1]] = ifelse(alpha[[1]]>=beta.ss, beta.ss-1, alpha[[1]])
beta[[1]] = beta.ss - alpha[[1]]

alpha[[2]] = p.a[[2]] + rnorm(R, 0, p.sigma)
alpha[[2]] = ifelse(alpha[[2]]>=beta.ss, beta.ss-1, alpha[[2]])
beta[[2]] = beta.ss - alpha[[2]]

#Lists to hold the results:
fluorescence = list()
fluorescence[[1]] = list()
fluorescence[[2]] = list()

p = list()
p[[1]] = list()
p[[2]] = list()

#Data-generating loop:
for (trt in 1:2) {
    for (reaction in 1:R) {
        #Initialize the reaction:
        X = rpois(1, lambda[[trt]])
        noise = rnorm(1,0,x.sigma)
        obs = ifelse(X<threshold, threshold+noise, X+noise)
        fluorescence[[trt]][[reaction]] = obs
        p[[trt]][[reaction]] = vector()     

        #Run the reaction chain:
        for (cycle in 1:N) {            
            p.k = rbeta(1, alpha[[trt]], beta[[trt]])
            X = X + rbinom(1, X, p.k)
            noise = max(threshold + rnorm(1,0,x.sigma), 0)
            obs = ifelse(X<threshold, threshold+noise, X+noise)
            fluorescence[[trt]][[reaction]] = c(fluorescence[[trt]][[reaction]], obs)
            p[[trt]][[reaction]] = c(p[[trt]][[reaction]], p.k)
            print(paste("cycle:", cycle))
        }
        print(paste("reaction:", reaction))
    }
    print(paste("trt:", trt))
}


