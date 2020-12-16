library(matrixStats)
library(hyper2)

run_bt <- function(df, player1="player1", player2="player2", iter=100000){
  player.names = sort(unique(c(df[,paste(player1)], df[,paste(player2)])))
  dataset.numerator = table(df$win)
  dataset.numerator = sapply(1:length(player.names), function(x){
    ifelse(player.names[x] %in% names(dataset.numerator), dataset.numerator[player.names[x]], 0)
  })
  names(dataset.numerator) = player.names
  #N <- length(player.names)
  N <- length(dataset.numerator)
  comparisons <- t(combn(N, 2))
  rownames(comparisons) <- 1:nrow(comparisons)
  #comparisons
  convert = data.frame(let = player.names, num = 1:length(player.names))
  comparisons_let = apply(apply(comparisons, c(1,2), function(x) paste(convert$let[match(x, convert$num)])), 1, 
                          function(y) paste(y, collapse=""))
  dataset.denominator = rep(0, length(comparisons_let))
  z = table(apply(df, 1, function(x) paste(sort(c(x[paste(player1)], x[paste(player2)])), collapse="")))
  dataset.denominator[match(names(z), comparisons_let)] = -z
  #encode the comparisons as a matrix
  comparisons.matrix <- sapply(1:N, function(x){apply(comparisons==x, 1, any)})
  comparisons.matrix <- as.matrix(sapply(1:N, function(x){t(matrix(apply(comparisons==x, 1, any)))}))
  
  if(nrow(comparisons.matrix) != nrow(comparisons)){
    comparisons.matrix = t(comparisons.matrix)
  }
  dataset <- hyper2(c(1:N,lapply(1:nrow(comparisons), function(x){
    comparisons[x,]})), c(dataset.numerator, dataset.denominator), pnames=player.names)
  
  ##########
  #MLE estimation using 'hyper2'
  #this gives the player ratings p, which are constrained to sum to 1
  mle.p <- maxp(dataset)
  mle.p
  
  ##########
  sample.crp.concentration <- function(){
    #sample latent variable eta
    zeta <- rbeta(1, alpha+1, N)
    ln.zeta <- log(zeta)
    #sample alpha conditional on eta
    z <- (prior.alpha.shape+k-1)/(N*(prior.alpha.rate-ln.zeta))
    pi <- z/(1+z)
    if (rbinom(1,1,pi)==1){
      alpha <- rgamma(1, prior.alpha.shape+k, prior.alpha.rate-ln.zeta)
    } else {
      alpha <- rgamma(1, prior.alpha.shape+k-1, prior.alpha.rate-ln.zeta)
    }
    alpha
  }
  m.rename <- function(m, string = T){
    unique.m <- unique(m)
    output <- sapply(m, function(y) {
      match(x = y, unique.m) - 1
    })
    if (string){
      paste(output, collapse = ",")
    } else {
      output
    }
  }
  ##########
  ##########
  #Gibbs sampler for Bayesian Bradley-Terry model
  #Caron and Doucet 2012, section 2: https://www.tandfonline.com/doi/full/10.1080/10618600.2012.638220
  #specify hyperparameter a; I think prior interpretation is p ~ Dirichlet(a)
  #not entirely sure what b does; I think something about the informativeness of the prior, but the paper says a=1 and b=0 leads to MAP=MLE
  a <- 1
  b <- 1
  #specify hyperparameters for prior concentration parameter
  #this assumes 50% probability that all are the same
  prior.alpha.shape <- 1
  prior.alpha.rate <- 2.333415
  #specify number of iterations for the sampler
  iter <- iter
  #create object to store results
  post.p <- matrix(NA, iter, N)
  post.C <- matrix(NA, iter, N)
  post.lambda <- matrix(NA, iter, N)
  post.alpha <- rep(NA, iter)
  colnames(post.p) <- player.names
  colnames(post.C) <- player.names
  colnames(post.lambda) <- player.names
  #set starting values
  C <- matrix(1, N, 1)
  C.list <- rep(1,N)
  p <- rep(1/N, N)
  theta <- a
  lambda <- theta[C.list]
  alpha <- prior.alpha.shape / prior.alpha.rate
  #precalculated quantities
  lgamma.a <- lgamma(a)
  #iterate sampler
  for (i in 1:iter){
    #print every 1000 iterations
    if (i%%1000==0){
      print(i)
    }
    #sample latent sum of lowest arrival times for each player comparison
    z <- rgamma(nrow(comparisons), -dataset.denominator, 
                sapply(1:nrow(comparisons), function(x){sum(lambda[comparisons[x,]], na.rm = T)})) #na.rm = T???
    #sample configuration of each player conditional on other players 
    for (j in 1:N){
      C[j,] <- 0
      lambda[j] <- NA
      if (!(1 %in% C[,C.list[j]])){
        C.update.index <- C.list[-j] > C.list[j]
        C.list[-j][C.update.index] <- C.list[-j][C.update.index] - 1
        C <- C[,-C.list[j],drop=F]
        theta <- theta[-C.list[j]]
      }
      C.list[j] <- NA
      k <- ncol(C)
      C.ln.ml <- sapply(1:k, function(x){dataset.numerator[j]*log(theta[x]) - theta[x]*sum(z[comparisons.matrix[,j]])})
      a.star <- a + dataset.numerator[j]
      #b.star <- 1 + sum(z[comparisons.matrix[,j]])
      b.star <- b + sum(z[comparisons.matrix[,j]])
      C.ln.ml <- c(C.ln.ml, -a.star*log(b.star) - lgamma.a + lgamma(a.star))
      C.ln.prior <- log(c(colSums2(C), alpha))
      C.prob <- C.ln.ml + C.ln.prior
      C.prob <- exp(C.prob - matrixStats::logSumExp(C.prob))
      C.indicator <- match(rmultinom(1,1,C.prob), x=1)
      C.list[j] <- C.indicator
      if (C.indicator > k){
        C <- cbind(C,0)
        C[j, C.indicator] <- 1
        k <- k + 1
        #sample new theta
        theta <- c(theta, rgamma(1, a.star, b.star))
      } else {
        C[j, C.indicator] <- 1
      }
      lambda[j] <- theta[C.list[j]]
    }
    #sample concentration parameter
    alpha <- sample.crp.concentration()
    #sample strengths for each class
    a.star <- a + dataset.numerator%*%C
    #b.star <- sapply(1:N, function(x){sum(z[comparisons.matrix[,x]])})%*%C
    b.star <- b + sapply(1:N, function(x){sum(z[comparisons.matrix[,x]])})%*%C
    theta <- rgamma(k, a.star, b.star)
    #resample total scale for strengths
    theta <- theta/sum(theta) * rgamma(1, k*a, 1)
    lambda <- theta[C.list]
    p <- lambda/sum(lambda)
    #impose order on groupings
    #C <- C[,order(-theta), drop=F]
    #C.list <- rank(-theta)[C.list]
    #theta <- theta[order(-theta)]
    #store results
    post.p[i,] <- p
    post.C[i,] <- m.rename(C.list, F)
    post.lambda[i,] <- lambda
    post.alpha[i] <- alpha
  }
  return(list(post.alpha=post.alpha, post.C=post.C, 
              post.p=post.p, post.lambda=post.lambda))
}

#compare marginal posterior means with MLE solution
#colMeans(post.p)
#mle.p

#calculate posterior probability that Karpov is the worst player (or tied for it)
#mean(apply(post.p, 1, which.min)==3)
#mean(apply(post.p, 1, min)==post.p[,3])
#sapply(1:ncol(post.p), function(x) mean(apply(post.p, 1, min)==post.p[,x]))

#quantiles for each player
#apply(post.p, 2, quantile, probs=c(0.025, 0.975))

#posterior number of groups
#table(apply(post.C,1,max))/iter

#posterior probability of each configuration
#head(-sort(-table(apply(post.C,1,paste,collapse=",")))/iter)

