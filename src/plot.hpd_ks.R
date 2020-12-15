source("graphics.R")

mcmc.stack <- function(coda.object){
  if(length(grep("array", class(coda.object))) > 0) {
    chains <- dim(coda.object)[3]
    all.mcmc <- sapply(1:dim(coda.object)[1], function(x){
          unlist(sapply(1:chains, function(y){
            as.vector(unlist(coda.object[x,,y]))}, simplify=F))
        }, simplify=T)
    colnames(all.mcmc) <- rownames(coda.object)
  }
  else {
    chains <- length(coda.object)
    all.mcmc <- list()
    if(chains == 1){
      all.mcmc[[1]] <- rbind(coda.object)
    } else {
      all.mcmc <- rbind(coda.object[[1]], coda.object[[2]])
      if (chains>2){
        for (j in 3:chains){
          all.mcmc <- rbind(all.mcmc, coda.object[[j]])
        }
      }
    }
    all.mcmc <- as.mcmc(all.mcmc)
  }
  
  return(all.mcmc)
}


get.inter.vars <- function(names)
{
  str <- matrix(NA, ncol=2, nrow=length(names))
  for(i in 1:length(names)){
    str[i,] <- strsplit(names[i], "[.]")[[1]]
  }
  return(str)
}


plot.reorder.ci <- function(med, mu,
                    hpd.narrow,
                    hpd.wide,
                    names=names,
                    xlab=xlab,
                    col.midvals="white",
                    pch.midvals="|",
                    type=type,
                    #pos=pos,
                    order=order,
                    remove=remove,
                    addline=addline,
                    ...)
{
  alldat <- data.frame(cbind(med, mu, hpd.narrow, hpd.wide))
  alldat$names <- factor(names, levels=names)
  removed <- gsub(paste0("^",remove),"",names)
  alldat$spec.names <- get.inter.vars(removed)
  if(order == 1){
    alldat$inter <- factor(alldat$spec.names[,1], levels=unique(alldat$spec.names[,1]))
  } else if (order == 2){
    alldat$inter <- factor(alldat$spec.names[,2], levels=unique(alldat$spec.names[,2]))
  } else {
    alldat$inter <- factor(alldat$names, levels=alldat$names)
  }
  ordered <- alldat[order(alldat$inter),]
  ordered$names <- factor(ordered$names, levels=rev(ordered$names))
  require(ggplot2)
  catplot <- ggplot(data=ordered) + 
          geom_point(aes(x=names, y=med), size=2) +
          geom_segment(aes(x=names, xend=names, y = lower.1, yend=upper.1), size=0.5) + 
          geom_segment(aes(x=names, xend=names, y = lower, yend=upper), size=1.5) +
          ylab("Effect estimate") + xlab("") +
          geom_point(aes(x=names, y=mu, group=inter), shape=pch.midvals, col=col.midvals, size=3) +
          coord_flip() + theme_minimal() 
  catplot <- catplot + geom_hline(yintercept = 0, colour="lightsteelblue", size=1, linetype="longdash") 
  if (addline) {
    catplot <- catplot + geom_line(aes(x=names, y=mu, group=inter), col="blue", size=1)
  } 
  
  return(catplot)          
}

plot.inter.ci <- function(med, mu,
                          hpd.narrow,
                          hpd.wide,
                          names,
                          xlab="HPD interval",
                          col.midvals="white",
                          pch.midvals="|",
                          grouped=rep(1,length(names)),
                          addline=F,
                          hpdgroup=rep(1,length(names)),
                          wide=T,
                          ordered=T, 
                          ylab = "Effect estimate",
                          title=NULL,
                          yint=0, 
                          ...)
{
  colnames(hpd.narrow) <- c("lower", "upper")
  colnames(hpd.wide) <- c("lower.1", "upper.1")
  
  tempdat <- data.frame(cbind(med, mu, hpd.narrow, hpd.wide), row.names = NULL)
  tempdat$names <- factor(rep(names, length.out=length(tempdat$mu)), levels=rev(names))
  tempdat$hpdgroup <- as.numeric(hpdgroup)
  tempdat$grouped <- as.character(rep(grouped, length.out=length(tempdat$mu)))
  
  if(ordered){
    alldat <- tempdat[order(tempdat$mu, decreasing = T), ]
    alldat$names <- factor(alldat$names, levels=unique(as.character(alldat$names)) )
  } else {
    alldat <- tempdat
  }

  require(ggplot2)
  alldat$colorLines <- factor(paste0(alldat$grouped,alldat$hpdgroup), levels=unique(paste0(alldat$grouped,alldat$hpdgroup)))
  catplot <- ggplot(data=alldat, aes()) + 
    geom_point(aes(x=names, y=mu, alpha=hpdgroup, col=grouped), size=3) + #col="#000000"
    geom_segment(aes(x=names, xend=names, y = lower, yend=upper, alpha=hpdgroup), col="#000000", size=1) + 
    ylab(ylab) + xlab("") + ggtitle(title) + 
    coord_flip() + theme_minimal() + 
    geom_hline(yintercept = yint, colour="gray50", size=0.5, linetype="longdash") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    # scale_colour_hue(guide=FALSE) + #palette="Blues", 
    scale_alpha_continuous(range=c(0.5,1), guide=FALSE)
  if (addline) {
    catplot <- catplot + geom_line(data=alldat, aes(x=names, y=mu, group=colorLines, 
                               colour=grouped, alpha=hpdgroup), size=0.5)
  } 
  if (wide) {
    catplot <- catplot + geom_segment(aes(x=names, xend=names, y = lower.1, yend=upper.1, alpha=hpdgroup), 
                                     col="#000000", size=0.5) 
  }
  
  catplot <- catplot + geom_point(aes(x=names, y=med), shape=pch.midvals, col=col.midvals, size=5)
  return(list(plot=catplot, df=alldat))         
}


plot.mine.hpd <- function(coda.object,
                     wanted=NULL,
                     prob.wide=0.95,
                     prob.narrow=0.50,
                     xlab="HPD interval",
                     names=NULL,
                     type="p",
                     centered=F,
                     title=NULL, 
                     #pos=0,  ## adjust the label's position
                     remove=NULL,
                     addline=F,
                     order=3,
                     grouped=rep(1,length(wanted)),
                     wide=T,
                     yint=0,
                     ...)
{
  chain <- list()
  mu <- c()
  med <- c()
  temphpd.wide <-list()
  temphpd.narrow <-list()
  if(is.null(wanted)) wanted=varnames(coda.object)
  if(is.null(varnames(coda.object))){
    if(length(grep("array", class(coda.object))) > 0){
      rownames(coda.object) = wanted
    } else {
      colnames(coda.object) = wanted
    }
  }
    

  #which.wanted=ifelse(is.integer(wanted), wanted, match(wanted, varnames(coda.object)))
  if((length(coda.object) > 1 && is.list(coda.object)) | length(grep("array", class(coda.object))) > 0) {
    chain <- as.mcmc(mcmc.stack(coda.object))
    nchain <- ifelse((length(coda.object) > 1 && is.list(coda.object)), length(chain), 1)
  } else {
    chain <- as.mcmc(coda.object)
    nchain <- ifelse(length(grep("list", class(coda.object))) > 0, length(chain), 1)
  }
  
  which.wanted=match(wanted, colnames(chain))
  num.wanted=length(which.wanted)
  for (i in 1:nchain){
    mu    <- c(mu, colMeans(chain[,which.wanted]))
    med   <- c(med, apply(chain[,which.wanted], 2, median))
    temphpd.wide[[i]]    <- coda::HPDinterval(chain, prob=prob.wide)[which.wanted,]
    temphpd.narrow[[i]]  <- coda::HPDinterval(chain, prob=prob.narrow)[which.wanted,]
  }
  #hpdgroup <- as.character(rep(1:length(chain), each=num.wanted))
  hpd.wide <- do.call(rbind, temphpd.wide)
  hpd.narrow <- do.call(rbind, temphpd.narrow)
  
  if (is.null(names)==T){
    names <- varnames(chain)[which.wanted]
  } else {
    names <- rep(names, length.out=num.wanted)
  }
  
  ypos <- plot.inter.ci(med, mu, hpd.narrow, hpd.wide, names=names, xlab=xlab, #hpdgroup=hpdgroup,
                        col.midvals="white", pch.midvals="|", addline=addline, grouped=grouped, 
                        wide=wide, title=title, yint=yint)#, ...)
                        #type=type, pos=pos, ...)

  invisible(ypos)
}


plot.hpd <- function(coda.object,
                     wanted=varnames(coda.object),
                     prob.wide=0.95,
                     prob.narrow=0.50,
                     xlab="HPD interval",
                     title=NULL, 
                     names=NULL,
                     type="p",
                     centered=F,
                     #pos=0,  ## adjust the label's position, 
                     cex.main=1,
                     name.lines = NULL,
                     ...)
{
  if (is.integer(wanted)){
    which.wanted <- wanted
    } else {which.wanted <- match(wanted, varnames(coda.object))}
  num.wanted=length(which.wanted)
  # chain <- mcmc.stack(coda.object)
  chain <- as.mcmc(coda.object)
  mu    <- colMeans(chain[,which.wanted])
  med   <- apply(coda::HPDinterval(chain, prob=0.01)[which.wanted,],
                 1, mean)
  hpd.wide    <- coda::HPDinterval(chain, prob=prob.wide)[which.wanted,]
  hpd.narrow  <- coda::HPDinterval(chain, prob=prob.narrow)[which.wanted,]
  
  if (is.null(names)) names <- varnames(chain)[which.wanted]
  else names <- rep(names, length.out=length(wanted))
  ypos <- plot.ci(med, hpd.narrow, hpd.wide, names=names, xlab=xlab, title=title, col.midvals="white", pch.midvals="|", cex.main=cex.main)
                  #,pos=pos, type=type, ...) +
  abline(v=0, col="gray")      
  if ("p"==type)
  {
    points(mu, ypos, pch="|")
  }
  invisible(ypos)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL, totTitle=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    title(totTitle)
    print(plots[[1]])
    title(totTitle, adj=1)
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

getPlot.compare.helper <- function(dfplot,ptypes,byVar)
{

  limz = c(-3,3)
  dfplot$pval[-log10(dfplot$pval)>(max(limz)-.01)] = 10^(-(max(limz)-.01))
  
  dfplot$dir.string = ""
  
  dfplot$p.with.direction = -log10(dfplot$pval) * dfplot$direction
  dfplot$stringg = 32
  dfplot$stringg[dfplot$pval<.05] = 8 ##symbol for sig.                                                                                                                                            
  
  aplot = ggplot(dfplot, aes(x=effectName1, y=effectName2, 
                             fill = p.with.direction, label=stringg, shape = stringg))
  aplot = aplot + geom_tile()
  aplot = aplot + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
  aplot = aplot + scale_shape_identity()
  aplot = aplot + xlab(paste(byVar,"1")) + ylab(paste(byVar,"2"))
  aplot = aplot + scale_fill_gradient2(low = "red", high = "blue",
                                       name = expression('-Log'[10]*'(p)'%.%'EffectDirection'),
                                       limits =limz)
  aplot = aplot + theme_classic() + geom_point(size = 2.5, color="white") + 
    facet_grid(Level~LevelPO) + ggtitle(paste(byVar,"contrast plot for",ptypes))
  aplot = aplot 
  return(aplot)
}

###############
# Ribbon plot #
###############

ribbon_plot <- function(mcmcList, ptypes, encoded, N=1, colors="rix", groups="poe", 
                        effects=c("rix","poe","diet"), grand_mean=0, 
                        Match=F, wide=0.95, narrow=0.5){
  savedata <- list()
  ribbon <- list()
  mu <- list()
  med <- list()
  #len <- 1
  its <- N
  mcmc_comb <- data.frame()
  for (j in 1:length(ptypes)){
    if(class(mcmcList) == "data.frame"){
      mcmcOb <- mcmcList
    } else if(length(mcmcList) == length(ptypes)){
      mcmcOb <- mcmcList[[j]]
    } else {
      n <- which(names(mcmcList) == ptypes[j])
      mcmcOb <- mcmcList[[n]]
    }
    savedata[[ptypes[j]]] <- list()
    for(k in 1:its){
      #if (class(mcmcOb) == "list"){
      #  its <- length(mcmcOb)
      each <- nrow(mcmcOb)/its
      low <- ((k-1)*each) + 1
      high <- low + each - 1
      mcmc_use <- as.mcmc(mcmcOb[low:high, ])
      #len <- length(mcmcOb)
      #}
      #else {
      #  mcmc_use <- as.mcmc(mcmcOb)
      #}
      mu[[k]]    <- c(colMeans(mcmc_use))
      med[[k]]   <- apply(coda::HPDinterval(mcmc_use, prob=0.01), 1, median)
      
      mcmc_comb <- rbind(mcmc_comb, mcmc_use)
    }
    
    mcmc_comb <- as.mcmc(mcmc_comb)
    mu_tot    <- c(colMeans(mcmc_comb))
    med_tot   <- apply(coda::HPDinterval(mcmc_comb, prob=0.01), 1, median)
    hpd.wide    <- coda::HPDinterval(mcmc_comb, prob=wide)
    hpd.narrow  <- coda::HPDinterval(mcmc_comb, prob=narrow)
    printdata <- data.frame(cbind(mu_tot, med_tot, hpd.wide,hpd.narrow))
    printdata$dietrix <- gsub("PO|_neg|_pos", "", rownames(printdata))
    #do.call(rbind,strsplit(rownames(printdata),"PO|_"))[,2] 
    printdata$diet <- gsub('[0-9]+', '', printdata$dietrix)
    printdata$rix <- gsub('[A-Z]+', '', toupper(printdata$dietrix))
    printdata <- printdata[which(printdata$dietrix %in% encoded$Level[which(encoded$Variable == "DietRIX")]),]
    if (Match){
      printdata <- printdata
    } else {
      printdata$poe <- "neg"
      printdata$poe[grep("PO|pos", rownames(printdata))] <- "pos"
      #do.call(rbind,strsplit(rownames(printdata),"PO"))[,1]
      printdata$poen <- ifelse(printdata$poe == "pos", 0.5, -0.5)
      printdata$rixpo <- paste0(printdata$rix, printdata$poe)
      #printdata <- printdata[-which(rownames(printdata) %in% c("STD10.neg","VDD10.neg","STD10.pos","VDD10.pos")),]
    }
    dietmeans <- aggregate(printdata[,1:2], list(printdata$diet), mean)
    dietmeans <- dietmeans[order(dietmeans[,"mu_tot"]),]
#browser()    
    printdata$diet <- factor(printdata$diet, levels = dietmeans$Group.1)
    printdata$rix <- factor(printdata$rix, levels = as.character(encoded$Level[which(encoded$Variable == "RIX")]))
    savedata[[ptypes[j]]][[k]] <- printdata
    ylab <- ifelse(length(grep("pos|neg", rownames(printdata)) > 0), "Predicted effect", "Estimated covariate effect")
    if (Match){
      ribbon[[ptypes[j]]] <- ggplot(printdata, aes(color=rix)) +
        geom_ribbon(color=NA, alpha=0.5, aes(x=diet, ymin=lower.1, ymax=upper.1,group=rix,
                                             fill=rix, color=rix)) + 
        geom_ribbon(color=NA, alpha=0.3, aes(x=diet, ymin=lower, ymax=upper,group=rix,
                                             fill=rix, color=rix)) + 
        geom_hline(yintercept = grand_mean, colour="lightsteelblue", size=0.5, linetype="longdash") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        scale_fill_discrete(guide=FALSE) + scale_color_discrete(guide=FALSE) + 
        facet_wrap(~rix) + theme_minimal() + 
        ylab(ylab) + xlab("") + ggtitle(ptypes[j])
      for(k in 1:its){
        #mu[[k]][match(names(mu[[k]]), rownames(printdata))]
        printdata$mu <- mu[[k]][match(names(mu[[k]]), rownames(printdata))][1:nrow(printdata)]
        ribbon[[ptypes[j]]] <- ribbon[[ptypes[j]]] + geom_line(data=printdata, aes(y=mu, x=diet, group=rix))
      }
      #ribbon[[ptypes[j]]] <- ribbon[[ptypes[j]]] + geom_line(data=printdata, aes(y=mu_tot, x=diet, group=rix), color="black")
    } else {
      ribbon[[ptypes[j]]] <- ggplot(printdata) + 
        geom_ribbon(colour=NA, aes_string(x=setdiff(effects, c(colors, groups)), ymin="lower", ymax="upper", 
                                          group=groups, alpha=groups, fill=colors, colour=colors)) +
        geom_ribbon(colour=NA, aes_string(x=setdiff(effects, c(colors, groups)), ymin="lower.1", ymax="upper.1", 
                                          group=groups, alpha=groups, fill=colors, colour=colors)) + 
        geom_line(aes_string(y="med_tot", x=setdiff(effects, c(colors, groups)), group=groups, alpha=groups), color="black") +  
        geom_hline(yintercept = 0, colour="lightsteelblue", size=0.5, linetype="longdash") +
        scale_alpha_discrete(range = rev(c(0.25, 0.65)), name="") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        scale_fill_discrete(guide=FALSE) + scale_color_discrete(guide=FALSE) + 
        facet_wrap(~rix) + theme_minimal() + 
        ylab(ylab) + xlab("") + ggtitle(ptypes[j])
    }
    print(paste("plots finished:", ptypes[j]))
  }
  return(list(ribbon=ribbon, savedata=savedata))
}



