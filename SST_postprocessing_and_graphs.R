# images and post processing for paper
rm(list=ls())

# Set-up directories
if(Sys.getenv('HOME')=='/home/scotti') {
  ## means it is on a server
  main_dir <- "/home/scotti/2022/behaveModelling" # LISA
  data_dir <- file.path(main_dir, 'SST/rdatas')
  res_dir <- file.path(main_dir, 'SST/samples')
  pdf_dir <- file.path(main_dir, 'SST/pdfs')
  cores <- 16
} else if (Sys.getenv('HOME')=='/Users/scotti') { # running locally using mountainDuck
  main_dir <- '/Users/scotti/Library/Group Containers/G69SCX94XU.duck/Library/Application\ Support/duck/Volumes/LISA/2022/behaveModelling'
  data_dir <- file.path(main_dir, 'SST/rdatas')
  res_dir <- file.path(main_dir, 'SST/samples')
  pdf_dir <- file.path(main_dir, 'SST/pdfs')  
  cores <- 4
}

setwd(main_dir)
source('dmc/dmc.R')
load_model ("EXG-SS","exgSSprobit.R")
load(file.path(res_dir, "img_save_ss.RData"))

##### HOW TO SAVE SUMMARY PARAMETER TABLE
params <- summary.dmc(hss)
tmp <- t(data.frame(lapply(params,function(x){x[[1]][,1]})))
tmp <- rbind(tmp,apply(tmp,2,mean))
row.names(tmp) <- c(names(hss),"Mean")
write.table(round(tmp,6), file='~/Downloads/dmc_params_final.tsv',sep='\t')

# set up parameter vector
p.vector  <- c(mu.true=.5,muS=.2,
               sigma.true=.05,sigmaS=.03,
               tau.true=.08,tauS=.05,
               tf=qnorm(.1),gf=qnorm(.1))

# Truncated normal priors:
p1 <- p.vector; p1[1:length(p1)] <- 1; p1[7:8] <- -1
p.prior <- prior.p.dmc(
  dists = rep("tnorm",length(p1)),p1=p.vector,p2=abs(p1), 
  lower=c(rep(0,length(p1)-2),-12,-12),upper=c(rep(2,length(p1)-2),rep(12,2)) # set uppers to 2 instead of 0.5 or 1
)
s.prior <- prior.p.dmc(p1=p1,p2=p1,
                       dists=rep("gamma",length(p1)))
pp.prior=list(p.prior,s.prior)
pp.prior <- attr(hss,"hyper")$pp.prior

# plot the priors
par(mfcol=c(2,4)); for (i in names(p.prior)) plot.prior(i,p.prior)

plot.dmc(hss,layout=c(3,6),p.prior=pp.prior,location=TRUE)
plot.dmc(hss,layout=c(3,6),p.prior=pp.prior,scale=TRUE)

# Create posterior predictives to check goodness of fit;
# This can be SLOW unless you have lots of cores available
pp <- h.post.predict.dmc(hss,cores=4, gglist=TRUE)

# We can now plot observed response proportions and the corresponding 
# 95% credible intervals of the posterior predictions 
# (i.e., the 2.5% and 97.5% quantiles of the predicted data).
# See lesson 5.5 for more detail on these functions. Because the proportion of
# non-responses on stop trials is of interest we must set the include.NR 
# argument to TRUE (by default proportions of trials with an RT response are
# plotted). The graph shows results at the group level, obtained by first
# generating predicted data for each subject and for each posterior sample, 
# then taking the average over subjects for each posterior sample. 
theme_set(theme_simple())
ggplot.RP.dmc(pp,include.NR=TRUE)  # Figure 13a

# Plot average cdf 
plot.pp.dmc(pp,layout=c(2,2),ylim=c(0,1),fits.pcol="grey",
            model.legend=FALSE,dname="",aname="",x.min.max=c(.25,2))

# The following plot functions require saving posterior predictve simulations;
# Again this can be SLOW
pp.sim <- h.post.predict.dmc(hss,save.simulation = TRUE,cores=4)

# To make these plots, we have to average over ranges of SSDs or else the result 
# is very noisy. This causes an issue when deriving averages because, as 
# illustrated above, different subjects can have very different SSD ranges.
# DMC averages in two ways, either in terms of absolute time or by dividing each
# subject's SSDs up according to percentile ranges. The following examples show
# that for the inhibition function the percentile range method is best (it is
# the default) as the absolute method flattens the function in a way that is 
# not representative of individual subjects. For the SRRT function, 
# the opposite it true, so the absolute method is the default.
par(mfrow=c(2,2))
plot_SS_if.dmc(data=hss,n.intervals=6,violin=.25,percentile.av=FALSE)
plot_SS_if.dmc(data=hss,n.intervals=6,violin=.25)
plot_SS_srrt.dmc(data=hss,n.intervals=6,violin=.25)
plot_SS_srrt.dmc(data=hss,n.intervals=6,violin=.25,percentile.av=TRUE)

# We can explicitly specify ranges with the probs argument (NB. you don’t need to specify 
# 0 or 1) or just provide n.intervals so you get evenly spaced ranges (here 
# 20%). Note that the posterior predictive p values can be sensitive to how intervals are defined.
plot_SS_if.dmc(data=hss,sim=pp.sim,n.intervals=5,violin=.25)

# Average absolute SSD for each percentile range
# (0,20%]  (20,40%]  (40,60%]  (60,80%] (80,100%] 
# 0.213     0.306     0.359     0.406     0.467 
# 
# Inhibition function 
# SSD   n    p
# 1   (0,20%] 295 0.01
# 2  (20,40%] 295 0.04
# 3  (40,60%] 295 0.91
# 4  (60,80%] 295 0.97
# 5 (80,100%] 295 0.56

# Here we use 9 percentile intervals. Note the 
# under-estimation at SSD = ~0.4s and over-estimation at SSD = ~.55s, confirming 
# the findings with the first set of intervals.
plot_SS_if.dmc(data=hss,sim=pp.sim,n.intervals=9,violin=.25,xlab="SSD")

# Average absolute SSD for each percentile range
# (0,11%]  (11,22%]  (22,33%]  (33,44%]  (44,56%]  (56,67%]  (67,78%]  (78,89%] (89,100%] 
# 0.186     0.259     0.300     0.331     0.360     0.386     0.408     0.440     0.483 
# 
# Inhibition function 
# SSD   n    p
# 1   (0,11%] 177 0.00
# 2  (11,22%] 148 0.02
# 3  (22,33%] 177 0.11
# 4  (33,44%] 147 0.59
# 5  (44,56%] 177 0.91
# 6  (56,67%] 148 0.63
# 7  (67,78%] 176 0.89
# 8  (78,89%] 148 0.97
# 9 (89,100%] 177 0.17

# Note the messy intervals due to random 
# perturbation of SSDs; for nice presentation, we can use xlabels to replace this.
# This is used to create Figure 14b
par(mfrow=c(1,1))
tmp <- plot_SS_srrt.dmc(data=hss,sim=pp.sim,n.intervals=9,
                        violin=.5,percentile.av=FALSE)
tmp <- round(attr(tmp,"cuts"),2)
xlabels <- paste("(",tmp[-length(tmp)],",",tmp[-1],"]",sep="")

plot_SS_srrt.dmc(data=hss,sim=pp.sim,n.intervals=9,
                 violin=.5,percentile.av=FALSE,xlabels=xlabels,xlab="SSD (s)")

# We can also plot correct RTs only:
par(mfrow=c(1,2))
plot_SS_srrt.dmc(data=hss,sim=pp.sim,n.intervals=9,violin=.25,
                 percentile.av=TRUE,do.correct=TRUE)
plot_SS_srrt.dmc(data=hss,sim=pp.sim,n.intervals=9,violin=.25,
                 do.correct=TRUE)

# There is little correlation among the population-level mean (location) parameters
location.r <- pairs.dmc(hss,location=TRUE,do.plot=FALSE)
summary(location.r)

# There is little correlation among the population-level standard deviation (scale) parameters
scale.r <- pairs.dmc(hss,scale=TRUE,do.plot=FALSE)
summary(scale.r)

### PRIORS FROM POSTERIORS ----

# The following function takes the posterior samples for one participant and
# uses them to update the p.prior object with parameters corresponding to the
# best fitting truncated normal distributions
post.prior1 <- make.tnorm.prior(p.prior,hss[[1]])

# Fitting is done by optimization and is not always guaranteed to work, so it
# is best to plot the results to check. In this case the fits are quite good
plot.dmc(hss[[1]],p.prior=post.prior1,layout=c(3,6))

# The following convenience function allows you to extract the numeric values 
# from the p.prior object as a named matrix
ppars <- get.ppars(post.prior1)
round(ppars,3)




# SAVE CDF, inhibition function and srrt image
pdf(file.path(pdf_dir,"cdf_if_srrt.pdf"),height=8,width = 9)
plot.pp.dmc(pp,layout=c(3,2),ylim=c(0,1),fits.pcol="grey",
            model.legend=FALSE,dname="",aname="",x.min.max=c(.25,2))
plot_SS_if.dmc(data=hss,sim=pp.sim,n.intervals=9,violin=.25,xlab="SSD")
plot_SS_srrt.dmc(data=hss,sim=pp.sim,n.intervals=9,
                 violin=.5,percentile.av=FALSE,xlabels=xlabels,xlab="SSD (s)")
dev.off()
