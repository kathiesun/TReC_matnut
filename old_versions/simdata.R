
library(data.table)

getTestCase <- function()
{

    n.diet = 4
    n.rix  = 9
    n.dam  = 200
    dam.pois = 6

    
    tau.rix          = 5
    tau.diet         = 5
    tau.rix.diet     = 1
    tau.rix.poe      = 1
    tau.rix.diet.poe = 1
    tau.dam          = 0 
    sigma            = 2

    ##TODO: incporporate cages, batches
    df.dam = data.table(data.frame(
        DamID  = factor(1:n.dam),
        Diet = factor(sample(1:n.diet, n.dam, replace = T)),
        RIX  = factor(sample(1:n.rix, n.dam, replace = T))),
        key = "DamID")

    
    df.dam[,PO:=sample(c(-.5, .5), .N, replace = T), by = c("RIX","Diet")]

    
    list(exists=.N>0)
    existing = df.dam[,.N,by=c("RIX","Diet", "PO")][, list(both = sum(PO>0) *sum(PO<0)), by = c("RIX","Diet")][both>0]
    existing$both = NULL
    df.dam = df.dam[existing, on=c("RIX","Diet")]
    df.dam[,n.pups := 1+rpois(1, lambda = dam.pois), by = "DamID" ]
    setkey(df.dam, "DamID")
    
    n = sum(df.dam$n.pups)
    df.pup  = data.table(data.frame(Pup.ID  = 1:n,
                                    epsilon = rnorm(n, mean = 0, sd = sigma),
                                    DamID     = factor(rep(df.dam$DamID, times = df.dam$n.pups))),
                         key = "DamID")

    #browser()
    

    df.dam$RIX.PO      = interaction(df.dam$RIX)
    df.dam$RIX.Diet    = interaction(df.dam$RIX, df.dam$Diet)
    df.dam$RIX.Diet.PO = interaction(df.dam$RIX, df.dam$Diet) 

    

    df.effects = rbind(
        data.table(variable = "tau.Diet",        level = NA, coef = tau.diet),
        data.table(variable = "tau.RIX" ,        level = NA, coef = tau.rix),
        data.table(variable = "tau.RIX.Diet" ,   level = NA, coef = tau.rix.diet),
        data.table(variable = "tau.RIX.PO" ,     level = NA, coef = tau.rix.poe),
        data.table(variable = "tau.RIX.Diet.PO", level = NA, coef = tau.rix.diet.poe),

        data.table(variable = "DamID",      level = 1:n.dam,  coef = rnorm(n.dam,   mean = 0, sd = tau.dam)), 
        data.table(variable = "Diet",       level = 1:n.diet, coef = rnorm(n.diet,  mean = 0, sd = tau.diet)),
        data.table(variable = "RIX",        level = 1:n.rix,  coef = rnorm(n.rix,   mean = 0, sd = tau.rix)),
        data.table(variable = "RIX.PO",     level = 1:n.rix,  coef = rnorm(n.rix,   mean = 0, sd = tau.rix.poe)),
        
        data.table(variable = "RIX.Diet",   level = do.call(paste, c(expand.grid(1:n.rix, 1:n.diet), sep = ".")),
                   coef = rnorm(n.rix*n.diet, mean = 0, sd = tau.rix.diet)),
        data.table(variable = "RIX.Diet.PO",level = do.call(paste, c(expand.grid(1:n.rix, 1:n.diet), sep = ".")),
                   coef = rnorm(n.rix*n.diet, mean = 0, sd = tau.rix.diet.poe)))

    df.dam$effect =
    df.effects[variable=="DamID"]       [df.dam,  on = c(level="DamID"       )]$coef +
    df.effects[variable=="Diet"]        [df.dam,  on = c(level="Diet"        )]$coef +
    df.effects[variable=="RIX"]         [df.dam,  on = c(level="RIX"         )]$coef +
    df.effects[variable=="RIX.PO"]      [df.dam,  on = c(level="RIX.PO"      )]$coef +
    df.effects[variable=="RIX.Diet"]    [df.dam,  on = c(level="RIX.Diet"    )]$coef + 
    df.effects[variable=="RIX.Diet.PO"] [df.dam,  on = c(level="RIX.Diet.PO" )]$coef

    df.pup = df.dam[df.pup]
    df.pup$outcome = df.pup$effect + df.pup$epsilon
    
    #df.pup$effect      = NULL
    df.pup$RIX.PO      = NULL
    df.pup$RIX.Diet    = NULL
    df.pup$RIX.Diet.PO = NULL
    return(list(cov.data = df.pup, effects = df.effects))
}

    
