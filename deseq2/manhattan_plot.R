library(lattice)

chr_lengths <- c(195323493, 130592705, 121970685, 120000968, 120318178, 124742797, 103899312,
                 98082450,  94869988, 90585865, 61305280, 182003861, 159916361, 156355905,
                 151702468, 149576561, 145324086, 129282844, 124493793, 170869960)
chr_lengths <- chr_lengths/1e6
manhattan.plot<-function(chr, pos, pvalue,
                         sig.level=NA, annotate=NULL, ann.default=list(), 
                         should.thin=T, thin.pos.places=2, thin.logp.places=2, 
                         main=NULL, xlab="Chromosome", ylab=expression(-log[10](p-value)),
                         region_list=NULL, 
                         col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {
  
  if (length(chr)==0) stop("chromosome vector is empty")
  if (length(pos)==0) stop("position vector is empty")
  if (length(pvalue)==0) stop("pvalue vector is empty")
  
  #make sure we have an ordered factor
  if(!is.ordered(chr)) {
    chr <- ordered(chr)
  } else {
    chr <- chr[,drop=T]
  }
  
  #make sure positions are in kbp
  if (any(pos>1e6)) pos<-pos/1e6;
  
  regions_list = lapply(regions_list, function(x){
    x$start_mb = x$start/1e6
    x$end_mb   = x$end/1e6
    return(x) 
  })
  
  #calculate absolute genomic position
  #from relative chromosomal positions
  posmin <- tapply(pos,chr, min);
  posmax <- tapply(pos,chr, max);
  posshift <- head(c(0,cumsum(posmax)),-1);
  names(posshift) <- levels(chr)
  genpos <- pos + posshift[chr];
  getGenPos<-function(cchr, cpos) {
    p<-posshift[as.character(cchr)]+cpos
    return(p)
  }
  regions_list = lapply(regions_list, function(x){
    x$start_gen = getGenPos(x$Chr, x$start_mb)
    x$end_gen = getGenPos(x$Chr, x$end_mb)
    return(x)
  })
  
  #parse annotations
  grp <- NULL
  ann.settings <- list()
  label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5, 
                      col=NULL, fontface=NULL, fontsize=NULL, show=F)
  parse.label<-function(rawval, groupname) {
    r<-list(text=groupname)
    if(is.logical(rawval)) {
      if(!rawval) {r$show <- F}
    } else if (is.character(rawval) || is.expression(rawval)) {
      if(nchar(rawval)>=1) {
        r$text <- rawval
      }
    } else if (is.list(rawval)) {
      r <- modifyList(r, rawval)
    }
    return(r)
  }
  
  if(!is.null(annotate)) {
    if (is.list(annotate)) {
      grp <- annotate[[1]]
    } else {
      grp <- annotate
    } 
    if (!is.factor(grp)) {
      grp <- factor(grp)
    }
  } else {
    grp <- factor(rep(1, times=length(pvalue)))
  }
#browser() 
  ann.settings<-vector("list", length(levels(grp)))
  ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)
  
  if (length(ann.settings)>1) { 
    lcols<-trellis.par.get("superpose.symbol")$col 
    lfills<-trellis.par.get("superpose.symbol")$fill
    for(i in 2:length(levels(grp))) {
      ann.settings[[i]]<-list(pch=pch, 
                              col=lcols[(i-2) %% length(lcols) +1 ], 
                              fill=lfills[(i-2) %% length(lfills) +1 ], 
                              cex=cex, label=label.default);
      ann.settings[[i]]$label$show <- T
    }
    names(ann.settings)<-levels(grp)
  }
  for(i in 1:length(ann.settings)) {
    if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
    ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
                                          parse.label(ann.settings[[i]]$label, levels(grp)[i]))
    #if(!labelPeaks & "label" %in% names(ann.settings[[i]])) ann.settings[[i]]$label$show = F
  }
  if(is.list(annotate) && length(annotate)>1) {
    user.cols <- 2:length(annotate)
    ann.cols <- c()
    if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
      ann.cols<-match(names(annotate)[-1], names(ann.settings))
    } else {
      ann.cols<-user.cols-1
    }
    for(i in seq_along(user.cols)) {
      if(!is.null(annotate[[user.cols[i]]]$label)) {
        annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, 
                                                    levels(grp)[ann.cols[i]])
      }
      ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], 
                                              annotate[[user.cols[i]]])
    }
  }
  rm(annotate)
  
  #reduce number of points plotted
  if(should.thin) {
    thinned <- unique(data.frame(
      logp=round(-log10(pvalue),thin.logp.places), 
      pos=round(genpos,thin.pos.places), 
      chr=chr,
      grp=grp)
    )
    logp <- thinned$logp
    genpos <- thinned$pos
    chr <- thinned$chr
    grp <- thinned$grp
    rm(thinned)
  } else {
    logp <- -log10(pvalue)
  }
  #rm(pos, pvalue)
  gc()
  
  #custom axis to print chromosome names
  axis.chr <- function(side,...) {
    if(side=="bottom") {
      panel.axis(side=side, outside=T,
                 at=((posmax+posmin)/2+posshift),
                 labels=levels(chr), 
                 ticks=F, rot=0,
                 check.overlap=F
      )
    } else if (side=="top" || side=="right") {
      panel.axis(side=side, draw.labels=F, ticks=F);
    }
    else {
      axis.default(side=side,...);
    }
  }
  lo = seq(floor(min(pos)), floor(min(pos))-10)[min(which(seq(floor(min(pos)), floor(min(pos))-10) %% 10 == 0))]
  hi = seq(ceiling(max(pos)), ceiling(max(pos))+10)[min(which(seq(ceiling(max(pos)), ceiling(max(pos))+10) %% 10 == 0))]
  
  
  axis.chr <- function(side,...) {
    if(side=="bottom") {
      panel.axis(side="bottom", outside=T,
                   at=seq(lo, hi, by=50),
                   labels=T,
                   ticks=T, rot=45,
                   check.overlap=T
      )
    } else if (side=="top" || side=="right") {
      panel.axis(side=side, draw.labels=F, ticks=F);
    } else {
      axis.default(side=side,...);
    }
  }
  
  #shade.blocks <- function(x, region_list){
  #  if(!is.null(region_list)){
  #    for(i in 1:length(region_list)){
  #      reg = region_list[[i]]
  #      for(j in 1:nrow(regions)){   
  #         panel.xblocks(x, x >= reg$Start[j],x <= reg$End[j], col = "lightgrey")
  #      }
  #    }
  #  }
  #}
  
  #make sure the y-lim covers the range (plus a bit more to look nice)
  prepanel.chr<-function(x,y,...) { 
    A<-list();
    maxy<-ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0)))+.5;
    A$ylim=c(0,maxy);
    A;
  }
  labs = names(ann.settings)[which(names(ann.settings) != "")]
  cols = unlist(lapply(labs, function(x) ann.settings[[x]]$col))
  if(length(cols) == 0) cols=seq(1:length(labs))
  xyplot(logp~pos | chr, groups=grp,#chr=chr, 
         axis=axis.chr, ann.settings=ann.settings, 
         prepanel=prepanel.chr, scales=list(relation="free", axs="i",x=list(at=pos)),
         key=list(columns=2, 
                  text=list(lab=labs),
                  points=list(pch=c(19,19), col=cols)),
         panel=function(x, y, ..., getgenpos) {
           if(!is.na(sig.level)) {
             #add significance line (if requested)
             panel.abline(h=-log10(sig.level), lty=2);
           }
           panel.superpose(x, y, ..., getgenpos=getgenpos);
           if(!is.null(panel.extra)) {
             panel.extra(x,y, getgenpos, ...)
           }
           if(!is.null(regions_list)){
             shade.blocks <- function(x=getgenpos, region_list){
               for(i in 1:length(region_list)){
                 reg = region_list[[i]]
                 for(j in 1:nrow(regions)){
                   panel.xblocks(x, x >= reg$Start[j],x <= reg$End[j], col = "lightgrey")
                 }
              }
             }
           }
         },
         panel.groups = function(x,y,..., subscripts, group.number) {
           A<-list(...)
           #allow for different annotation settings
           gs <- ann.settings[[group.number]]
           A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
           A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
           A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
           A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
           A$x <- x
           A$y <- y
           do.call("panel.xyplot", A)
           #draw labels (if requested)
           if(gs$label$show) {
             gt<-gs$label
             names(gt)[which(names(gt)=="text")]<-"labels"
             gt$show<-NULL
             if(is.character(gt$x) | is.character(gt$y)) {
               peak = which.max(y)
               center = mean(range(x))
               if (is.character(gt$x)) {
                 if(gt$x=="peak") {gt$x<-x[peak]}
                 if(gt$x=="center") {gt$x<-center}
               }
               if (is.character(gt$y)) {
                 if(gt$y=="peak") {gt$y<-y[peak]}
               }
             }
             if(is.list(gt$x)) {
               gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
             }
             do.call("panel.text", gt)
           }
         },
         xlab=xlab, ylab=ylab, main=main, 
         panel.extra=panel.extra, getgenpos=getGenPos, ...
  );
}


genome_wide_plot <- function(chr, pos, pvalue,
                             sig.level=NA, annotate=NULL, regions=NULL, chr_lengths=NULL, 
                             main=NULL, 
                             xlab="Chromosome", ylab=expression(-log[10](p-value)),
                             col=c("seagreen1")){
  rix_mat = data.frame(pos = pos/1e6,
                       chr = chr,
                       pvalue = -log10(pvalue),
                       ann = ann)
  #rix_mat$chr <- unlist(lapply(rix_mat$chromosome, function(x) ifelse(x == "X", 20, as.numeric(x))))
  #if(!is.null(chr_lengths)) rix_mat$add_pos = chr_lengths[1:nrow(rix_mat)]
  rix_mat %>% arrange(chr, pos) %>% 
    group_by(chr) %>% dplyr::slice(n()) %>% 
    dplyr::select(one_of("chr", "pos")) %>%
    group_by() %>%
    #mutate(add_pos = pos, add_pos = if(!is.null(chr_lengths)) chr_lengths[1:length(unique(rix_mat$chr))]) %>%
    mutate(cum_pos = cumsum(pos)) %>% 
    group_by() %>%
    arrange(chr) -> chr_pos
  #if(!is.null(chr_lengths)) chr_pos$cum_pos
  
  half_pos <- unlist(sapply(1:nrow(chr_pos), function(i){
    tmp = ifelse(length(chr_pos$cum_pos[i-1]) > 0, chr_pos$cum_pos[i-1], 0)
    (chr_pos$cum_pos[i] + tmp)/2
  }))

  regions$new_start = apply(regions, 1, function(x) 
    suppressWarnings(
      ifelse(as.numeric(x["Chr"]) > 1, as.numeric(x["start"])/1e6 + chr_pos$cum_pos[which(chr_pos$chr %in% x["Chr"])-1], 
             as.numeric(x["start"])/1e6))
  )
  regions$new_end = apply(regions, 1, function(x) 
    suppressWarnings(
      ifelse(as.numeric(x["Chr"]) > 1, as.numeric(x["end"])/1e6 + chr_pos$cum_pos[which(chr_pos$chr %in% x["Chr"])-1], 
             as.numeric(x["end"])/1e6))
  )
  rix_mat %>%
    group_by(chr) %>%
    arrange(chr, pos) -> rix_mat
  new_pos = apply(rix_mat, 1, function(x) 
    ifelse(as.numeric(x["chr"]) > 1, as.numeric(x["pos"]) + chr_pos$cum_pos[which(chr_pos$chr == as.numeric(x["chr"]))-1], x["pos"])
  )
  rix_mat$plot_pos = as.numeric(paste(new_pos))
  rix_mat <- data.frame(rix_mat) %>% arrange(ann)
  
  y_max = max(rix_mat$pvalue, na.rm=T)
  xranges = seq(min(rix_mat$pos), max(rix_mat$pos))
  plot(x=xranges, y=rep(0, length(xranges)), main=main, 
       xlab=xlab, ylab=ylab, lty="blank",
       xlim=c(0,max(rix_mat$plot_pos)),ylim=c(0, y_max+y_max*0.1), 
       xaxs="i", yaxs="i", xaxt="n")
  axis(1, at=half_pos, labels=unique(rix_mat$chr))
  for(c in unique(rix_mat$chr)){
    c = as.numeric(c)
    shade_col = ifelse(c%%2 == 0, "gray85", "white")
    rix_mat %>% filter(chr == c) %>% 
      arrange(pos)-> ranges
    polygon(x=c(ranges$plot_pos[1], ranges$plot_pos[1], ranges$plot_pos[nrow(ranges)],ranges$plot_pos[nrow(ranges)]),
            y=c(-1,y_max+y_max*0.5,y_max+y_max*0.5,-1), 
            border = shade_col,
            col=shade_col) 
  }
  for(seg in 1:nrow(regions)){
    polygon(x=c(regions[seg, "new_start"], regions[seg, "new_start"], 
                regions[seg, "new_end"], regions[seg, "new_end"]),
            y=c(-1,y_max+y_max*0.5,y_max+y_max*0.5,-1), 
            col=rgb(0.12,0.56,1.00, 0.2), border=NA)
  }
  points(x=rix_mat$plot_pos, y=rix_mat$pvalue, 
         col=ifelse(rix_mat$ann == "", "slategray4", col[1]), pch=20)
  
}



