M = 10^4 # number of Monte Carlo simulations used to approximate integrals.

func1 = function(z_0){ # this function simulates from the distribution of the first group N(-2, 1/3) and returns the average probability of connection from nodes in this group to the reference node at (0, 0)

  z = matrix(rnorm(2 * M, 2, 1/3), ncol = 2)

  dist.temp = apply(t(z), 2, function(center) {

    norm(z_0 - center, "2")

  })

  values = exp(-dist.temp)

  return(mean(values))

}



func2 = function(z_0){ # this function simulates from the distribution of the second group N(-2, 1/3) and returns the average probability of connection from nodes in this group to the reference node at (0, 0)

  z = matrix(rnorm(2 * M, -2, 1/3), ncol = 2)

  dist.temp = apply(t(z), 2, function(center) {

    norm(z_0 - center, "2")

  })

  

  values = exp(-dist.temp)

  return(mean(values))

}





func = function(z_0){ # this function returns the difference between the ratio of link probabilities from the reference node to groups 1 and 2 and substracts the observed frequency of edges to groups 1 and 2. This is the objective function whose zero is the estimate of the node location. 

  return(abs(func1(z_0) / func2(z_0) - exp(a)))

}





z_0 = c(0, 0) 

numSim = 25

n.vec = c(50, 100, 500, 1000, 10^4)



diff = matrix(rep(0, length(n.vec)  * numSim), ncol = length(n.vec))

z.hat.saved = array(rep(0, 2 * numSim * length(n.vec) ), dim = c(2, numSim, length(n.vec)))



for(indexSim in 1:numSim){

  for(indexN in 1:length(n.vec)){

    print(print(indexSim))

    n = n.vec[indexN]

    group_1_location = matrix(rnorm(2 * n, 2, 1/3), ncol = 2)

    group_2_location = matrix(rnorm(2 * n, -2, 1/3), ncol = 2)

    

    dist_1 = apply(t(group_1_location), 2, function(point) {

      norm(z_0 - point, "2")

    })

    

    dist_2 = apply(t(group_2_location), 2, function(point) {

      norm(z_0 - point, "2")

    })

    

    

    P_group_1 = exp(-dist_1)

    P_group_2 = exp(-dist_2)

    

    Y_i1 = Y_i2 = 0

    while(sum(Y_i1) == 0 || sum(Y_i2) == 0){

      Y_i1 = rbinom(n, 1, P_group_1)

      Y_i2 = rbinom(n, 1, P_group_2)

    }

    

    a = log(mean(Y_i1)) - log(mean(Y_i2))

    z.hat = optim(c(0, 0), func, method = "L-BFGS-B")$par

    z.hat.saved[,indexSim, indexN] = z.hat

    diff[indexSim,indexN] = norm(z.hat, "2")

  }

}





pdf("compare_n_z_i_scatter.pdf")

par(mfrow=c(1,2))    # set the plotting area into a 1*2 array



plot(z.hat.saved[1,,1], z.hat.saved[2,,1], pch = 19, xlab = "", 

     ylab = "", main=TeX('Estimate of node location $\\z_i$ with n = 50'))

plot(z.hat.saved[1,,5], z.hat.saved[2,,5], pch = 19,  xlab = "", 

     ylab = "", main=TeX('Estimate of node location $\\z_i$ with n = 10^4'))

dev.off()











# Now let's estimate the node effects 

estimate_nu = function(n, numSim, nu.0){

  nu.hat = rep(0, numSim)

  nu = runif(n, -2, 0)

  # Simulate the expressions E[exp{-d(z, z_i)}] and E(exp(nu))

  group.location.temp = matrix(rnorm(2 * 10^4, 2, 1/3), ncol = 2) # create a temporary set of locations

  dist.temp = apply(t(group.location.temp), 2, function(point) {norm(z_0 - point, "2")})  # compute the distance between z_0 and group.location.temp

  tau = log(mean(exp(-dist.temp))) # estimate log(E(exp(-d(z, z_i))))                             

  group.1.location = matrix(rnorm(2 * n, 2, 1/3), ncol = 2) # simulate locations of the nodes in group k

  dist = apply(t(group.1.location), 2, function(point) { norm(z_0 - point, "2") } ) # compute distance between z_0 and group.1.location

  P = exp(-dist + nu.0 + nu)       # compute probability of edge

  

  for(indexSim in 1:numSim){

    Y_i1 = sum(rbinom(n, 1, P))          # simulate ARD

    nu.hat[indexSim] = log(Y_i1/n) - tau - log(mean(exp(nu))) # compute estimate of v_i

  }

  return(nu.hat)

}



n.vec = c(250, 500, 1000, 10^4)

numSim = 100

nu.0 = -1

nu.estimates = matrix(rep(0, length(n.vec) * numSim), ncol = length(n.vec))

for(index.n in 1:length(n.vec)){

  nu.estimates[,index.n] = estimate_nu(n.vec[index.n], numSim, nu.0)

}



pdf("node_effects.pdf")

boxplot(nu.estimates, names = n.vec, cex.lab = 1.5, cex.axis = 1.25, 

        xlab = "n", ylab = "Estimate of node effect", main=TeX('Estimate of node effect $\\nu_i$'))

abline(h =-1, col="red", lwd=3, lty=2)

dev.off()
