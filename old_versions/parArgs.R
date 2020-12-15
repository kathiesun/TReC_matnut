##Load files up (fread reccomended), get them into the right form.
#Note that doreper's fitBoxCox expects each column of y.mat to be a different phenotype, NOT each row.

##~30K phenotypes, so probably want to do this on the cluster; set up the parallel args for
##fit box cox
getBestParArgs <- function(batchSize=100, mem.gb = 4, timeLimit.hours = 23, onCluster=T, system.type="killdevil", filesToSource)
{
  if(onCluster)
  {
    parallelArgs = list(system.type     = system.type,
                        filesToSource   = filesToSource,
                        batchSize       = batchSize,
                        timeLimit.hours = timeLimit.hours,
                        cpuMemLimit.GB  = mem.gb,
                        outdir          = "/netscr/kys6/tmp",
                        saveProp        = T,
                        coresPerJob     = 1)
    
  } else {
    parallelArgs = getLocalParArgs(1)
  }
  
  return(parallelArgs)
}

getLocalParArgs <- function(mc.cores)
{
    parallelArgs = list(mc.cores = mc.cores,
                        mclBatch = mc.cores*100)
    return(parallelArgs)
    
}

