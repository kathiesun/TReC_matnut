options(digits=4)

setwd("C:/Users/Kathie/TReC_matnut/src")
library(tidyverse)
library(qtl2)

dir <- "C:/Users/Kathie/Dropbox\ (ValdarLab)"

##################################
##    heterozygous ideograms    ##
##################################

## allele probabilities for mnt and cegs
alleleprobs_minimuga <- readRDS(file.path(dir, "mini/interp_rqtl_allChr_alleleprob_17mar2020.rds"))
alleleprobs_4 = alleleprobs_minimuga[["4"]]
alleleprobs_7 = alleleprobs_minimuga[["12"]]
alleleprobs_12 = alleleprobs_minimuga[["12"]]
alleleprobs_minimuga = NULL
## all markers with probabilities
markers_7 = dimnames(alleleprobs_7)[[3]]

all_genes = read.csv(file.path(dir, "mini/combined_map_kmerMarker_17mar2020.csv"))

## make probability objects into dataframes and combine
alleleprobs_7_df = list()
for(i in 1:dim(alleleprobs_7)[1]){
  tmp = t(alleleprobs_7[i,,])
  tmp = data.frame(marker = rownames(tmp), tmp) %>% left_join(all_genes, by="marker")
  alleleprobs_7_df[[dimnames(alleleprobs_7)[[1]][i]]] = tmp
}


alleleprobs_7_df = do.call("rbind", alleleprobs_7_df)
alleleprobs_7_df$Pup.ID = as.integer(unlist(strsplit(rownames(alleleprobs_7_df),"[.]"))[c(T,F)])
alleleprobs_7_df = alleleprobs_7_df %>% 
  left_join(lab, by="Pup.ID")

## sort through dataframe and remove missing data
rem = which(is.nan(alleleprobs_7_df$A))
if(length(rem) > 0) { 
  alleleprobs_7_df[-rem, ] %>% arrange(RRIX, Pup.ID, pos) %>%
    distinct() -> alleleprobs_7_df
}

alleleprobs_7_df %>%
  filter(!Pup.ID %in% problemPups) %>%
  arrange(RRIX, pos, Pup.ID) %>% distinct() -> alleleprobs_7_df

## round out probs
alleleprobs_7_df[,LETTERS[1:8]] = round(alleleprobs_7_df[,LETTERS[1:8]], 3)
alleleprobs_7_df$RIX = as.factor(paste0(alleleprobs_7_df$CC.1,"/",alleleprobs_7_df$CC.2))
alleleprobs_7_df = alleleprobs_7_df %>% arrange(RRIX, pos) %>% select(-X)
keep_round = alleleprobs_7_df
keep_round[,LETTERS[1:8]] = round(keep_round[,LETTERS[1:8]], 1)
keep_round = keep_round %>% select(-"Pup.ID") %>% 
  arrange(RIX, pos) %>% distinct() %>%
  mutate(RIX_pos = paste0(pos, "_", RIX))

## only keep positions with agreement among all samples in a CC-RIX
#keep_pos = names(table(keep_round$RIX_pos))[which(table(keep_round$RIX_pos) == 1)]
keep_pos = keep_round$RIX_pos
alleleprobs_comb_nod = alleleprobs_7_df %>% select(-"Pup.ID") %>% 
  arrange(RIX, pos) %>% distinct() %>%
  mutate(RIX_pos = paste0(pos, "_", RIX)) %>%
  filter(RIX_pos %in% keep_pos) %>%
  group_by(RRIX, RIX, RIX_pos, CC.1, CC.2, marker, pos) %>%
  dplyr::summarize(A = mean(A),B = mean(B),C = mean(C),D = mean(D),
                   E = mean(E),F = mean(F),G = mean(G),H = mean(H))

## identify homozygous regions 
homog = apply(alleleprobs_comb_nod, 1, function(x) ifelse(any(x[LETTERS[1:8]] > 0.8), T, F))
plot_hets=data.frame(alleleprobs_comb_nod)
plot_hets$homog = homog

## determine top two founders
het_founds = data.frame(do.call("rbind",apply(plot_hets, 1, function(x){
  sort(names(sort(x[LETTERS[1:8]][which(x[LETTERS[1:8]] > 0.05)],decreasing=T)[1:2]))
})))

colnames(het_founds) = c("founder1", "founder2")
plot_hets = cbind(plot_hets, het_founds)
plot_hets$founders = apply(plot_hets, 1, function(x) 
  paste(sort(c(paste(x["founder1"]), paste(x["founder2"]))), collapse = ""))
plot_hets %>% 
  arrange(RIX, pos) %>%
  group_by(founders) -> plot_hets

## wide-to-long format
plot_hets_points = plot_hets %>% select("RIX","pos","founder1","founder2","homog") %>%  
  gather(which_f, founder, founder1:founder2, factor_key = T) %>% 
  #rename("Founder"="founder") %>%
  arrange(RIX, pos, which_f) 

## order in plot

plot_hets_points$RIX = factor(plot_hets_points$RIX, levels=unique(plot_hets_points$RIX), ordered=T)

## convert x-axis to numeric such that each CC-RIX gets two points
plot_hets_points$Pup = as.numeric(factor(plot_hets_points$RIX)) +  
  as.numeric(paste0("0.",gsub("founder","", as.character(plot_hets_points$which_f))))*2
#plot_hets_points$Founder = factor(plot_hets_points$Founder, levels = LETTERS[1:8])
npup = length(unique(plot_hets_points$RIX))


## lower transparency of homozygous regions
plot_hets_points$alpha = 1
plot_hets_points$alpha[which(plot_hets_points$homog)] = 0.0001



low = 59
hi1=hi2= 63
ribbon = data.frame(x=c(-Inf, Inf),
                    lo=low,  hi=min(hi1, hi2))

## classic CC colors
cols = c("yellow", "gray", "pink", "blue", "deepskyblue", "forestgreen", "red", "purple")

p = ggplot() + 
  geom_ribbon(data=ribbon, aes(ymin=lo, ymax=hi, x=x), fill="#a3752e", alpha=0.2) + 
  geom_point(data=plot_hets_points, aes(x=Pup, y=pos, col=founder, alpha=alpha), size=0.5,shape=15) + 
  scale_colour_manual(values=cols) + 
  scale_alpha_continuous(NULL, NULL) + 
  scale_x_continuous(name ="CC-RIX", breaks=c(1:npup)+0.3, #96
                     labels = as.character(unique(plot_hets_points$RIX)),#Pup.ID
                     limits=c(1,npup+0.5)) +
  scale_y_continuous(name ="Position", breaks=seq(0,170,5)) +
  theme_classic() + 
  theme(axis.text.x=element_text(angle=45,hjust=1, size=8),
        axis.text.y=element_text(size=8))

p
