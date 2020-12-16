#package for maximum likelihood estimation of Bradley-Terry model
library(hyper2)

##########
#specifying example data 'chess' from 'hyper2'; vignette table 1
#'hyper2' vignette: https://cran.r-project.org/web/packages/hyper2/vignettes/hyper2.pdf
#specify player names
player.names <- c("Topalov", "Anand", "Karpov")
#numerator for the likelihood; TOTAL number of times each player won versus ANY player
dataset.numerator <- c(30, 36, 22)
#enumerate all possible player comparisons
N <- length(dataset.numerator)
comparisons <- t(combn(N, 2))
rownames(comparisons) <- 1:nrow(comparisons)
comparisons
#encode the comparisons as a matrix
comparisons.matrix <- sapply(1:N, function(x){apply(comparisons==x, 1, any)})
comparisons.matrix
#denominator for the likelihood; NEGATIVE number of times each comparison occurred
dataset.denominator <- -c(35, 18, 35)
##########

#frequentist implementation of Bradley-Terry model
#format data as 'hyper2' object
dataset <- hyper2(c(1:N,lapply(1:nrow(comparisons), function(x){comparisons[x,]})), c(dataset.numerator, dataset.denominator), pnames=player.names)
dataset
#MLE estimation using 'hyper2'
#this gives the player ratings p, which are constrained to sum to 1
mle.p <- maxp(dataset)
mle.p

##########
#Gibbs sampler for Bayesian Bradley-Terry model
#Caron and Doucet 2012, section 2: https://www.tandfonline.com/doi/full/10.1080/10618600.2012.638220
#specify hyperparameter a; I think prior interpretation is p ~ Dirichlet(a)
#not entirely sure what b does, I think something about the informativeness of the prior, but the paper says a=1 and b=0 leads to MAP=MLE
a <- 1
b <- 0
#specify number of iterations for the sampler
iter <- 100000
#create object to store results
post.p <- matrix(NA, iter, N)
colnames(post.p) <- player.names
#set starting values based on MLE solution
p <- mle.p
names(p) <- NULL
lambda <- a*N*p
#store posterior hyperparameter a.star
a.star <- a + dataset.numerator
#iterate sampler
for (i in 1:iter){
  #print every 1000 iterations
  if (i%%1000==0){
    print(i)
  }
  #sample latent variable z
  z <- rgamma(nrow(comparisons), -dataset.denominator, sapply(1:nrow(comparisons), function(x){sum(lambda[comparisons[x,]])}))
  #update posterior hyperparameter b.star and unconstrained player ratings lambda
  b.star <- b + sapply(1:N, function(x){sum(z[comparisons.matrix[,x]])})
  lambda <- rgamma(N, a.star, b.star)
  #calculate p given lambda
  p <- lambda/sum(lambda)
  #resample lambda given p
  lambda <- p * rgamma(1, N*a, 1)
  #store posterior sample
  post.p[i,] <- p
}
#compare marginal posterior means with MLE solution
colMeans(post.p)
mle.p
#calculate posterior probability that Karpov is the worst player
mean(apply(post.p, 1, which.min)==3)
#quantiles for each player
apply(post.p, 2, quantile, probs=c(0.025, 0.975))