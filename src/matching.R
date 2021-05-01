matching = new.env()

##matchon: a field we want to be identical between pairs
matching$generateMatching <- function(df, matchon, matchoff, idcol)
{
    dfbad = list()
    df = data.table(df)
##    dfnew = df[,j=list(ID = paste(get(idcol), collapse=",")), by = c(matchon, matchoff)]
    dfnew = df[,j=list(ID = paste(get(idcol), collapse=",")), by = c(matchon)]

    l1.ID = c()
    l2.ID = c()
    offlevels = unique(df[[matchoff]])
    for(i in 1:nrow(dfnew))
    {
        row = dfnew[i,]
        IDs = unlist(strsplit(row$ID, ","))
        matchingOffValues = df[[matchoff]][match(IDs, df[[idcol]])]

        level1.IDs  = IDs[which(matchingOffValues == offlevels[1])]
        level2.IDs  = IDs[which(matchingOffValues == offlevels[2])]

        
        numpairs = min(length(level1.IDs), length(level2.IDs))
        if(numpairs == 0)
        {
            dfbad = util$appendToList(dfbad, dfnew[i,])
            next
        }
        
        level1.IDs = sample(x = level1.IDs, replace = F,size = numpairs)
        level2.IDs = sample(x = level2.IDs, replace = F,size = numpairs)

        l1.ID = c(l1.ID, level1.IDs)
        l2.ID = c(l2.ID, level2.IDs)
    }
    out = data.frame(ID.1 = l1.ID, ID.2 = l2.ID)
    for(cname in matchon)
    {
        for(j in 1:2)
        {
            out[[paste0(cname,".",j)]] = df[[cname]][match(l1.ID, df[[idcol]])]
        }
    }

    dfbad = do.call(rbind, dfbad)
    return(list(out=out, dfbad=dfbad))
}

##matchon: a field we want to be identical between pairs
matching$generateEveryMatch <- function(df, matchon, matchoff, idcol)
{
    dfbad = list()
    df = data.table(df)
    ##    dfnew = df[,j=list(ID = paste(get(idcol), collapse=",")), by = c(matchon, matchoff)]
    dfnew = df[,j=list(ID = paste(get(idcol), collapse=",")), by = c(matchon)]
    
    l1.ID = c()
    l2.ID = c()
    offlevels = unique(df[[matchoff]])
    for(i in 1:nrow(dfnew))
    {
        row = dfnew[i,]
        IDs = unlist(strsplit(row$ID, ","))
        matchingOffValues = df[[matchoff]][match(IDs, df[[idcol]])]
        
        level1.IDs  = IDs[which(matchingOffValues == offlevels[1])]
        level2.IDs  = IDs[which(matchingOffValues == offlevels[2])]
        
        
        numpairs = prod(length(level1.IDs), length(level2.IDs))
        if(numpairs == 0)
        {
            dfbad = util$appendToList(dfbad, dfnew[i,])
            next
        }
        combs = data.frame(expand.grid(level1.IDs, level2.IDs))
        level1.IDs = paste(combs$Var1)
        level2.IDs = paste(combs$Var2)
        
        l1.ID = c(l1.ID, level1.IDs)
        l2.ID = c(l2.ID, level2.IDs)
    }
    out = data.frame(ID.1 = l1.ID, ID.2 = l2.ID)
    for(cname in matchon)
    {
        for(j in 1:2)
        {
            out[[paste0(cname,".",j)]] = df[[cname]][match(l1.ID, df[[idcol]])]
        }
    }
    
    dfbad = do.call(rbind, dfbad)
    return(list(out=out, dfbad=dfbad))
}




matching$getDeltaForPhen <- function(original, matched, idcol, phen)
{
    y1 = original[[phen]][match(matched$ID.1, original[[idcol]])]
    y2 = original[[phen]][match(matched$ID.2, original[[idcol]])]
    out = y2 - y1
    matched$y = out
    return(matched)
}
