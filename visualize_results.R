# script to analyse results of SST modelling
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

# set some variables
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

# look at mean summary measures per person per parameter
summary.dmc(hss)

# plot posterior likelihoods
plot.dmc(hss,hyper=TRUE,pll.chain=TRUE)

# Check chains for the population-level (hyper) parameters
plot.dmc(hss,hyper=TRUE,layout=c(3,6))

# Get info criterion
h.IC.dmc(hss,DIC=TRUE)

# check gelman diags (results <1.1 indicate convergence)
# evaluates MCMC convergence by analyzing the difference between multiple Markov chains
# The convergence is assessed by comparing the estimated between-chains and within-chain variances for each model parameter. 
# Large differences between these variances indicate nonconvergence
gelman.diag.dmc(hss)
# check hyper as well
gelman.diag.dmc(hss,hyper=TRUE) 

# Plot CDFs for left/rihgt on go/ss trials
plot.pp.dmc(ppss,model.legend = TRUE,layout=c(2,2))

# pdf(file.path(pdf_dir,'SST_model_fit.pdf'),width = 13,height = 9) # make PDF file of graphs
# Plot CDFs for left/rihgt on go/ss trials
plot.pp.dmc(ppss,model.legend = TRUE,layout=c(2,2), dname='Data', aname='',style='cdf')
# dev.off()
# plot median signal-respond RTs
plot_SS_srrt.dmc(hss)

# plot overall inhibbtion function
plot_SS_if.dmc(hss, n.intervals=7,violin=.25,sim=ppss)

# plot hyper prior distributions for the popultion means
par(mfcol=c(2,4)); for (i in names(p.prior)) plot.prior(i,p.prior)

# plot hyper propr and posterior dists for the pop means
plot.dmc(hss,layout=c(3,6),p.prior=pp.prior,density = TRUE,prior.col = 'red')
plot.dmc(hss,hyper=TRUE,layout=c(3,3),p.prior=pp.prior)

# plot Hyper-prior and posterior distributions for the population standard deviations
plot.dmc(hss,layout=c(3,6),p.prior=pp.prior,scale=TRUE)


#####
#######
# IMAGES TO SAVE
########

pdf(file.path(pdf_dir,'SST_model_fit.pdf'),width = 13,height = 9) # make PDF file of graphs
# Plot CDFs for left/rihgt on go/ss trials
plot.pp.dmc(ppss,model.legend = TRUE,layout=c(2,2), dname='Data', aname='',style='cdf')
dev.off()

pdf(file.path(pdf_dir,'SST_model_priors.pdf'),width = 6,height = 5.5) # make PDF file of graphs
par(mfcol=c(3,3)); for (i in names(p.prior)) plot.prior(i,p.prior)
dev.off()

 pdf(file.path(pdf_dir,'SST_model_chains.pdf'),width = 10,height = 6) # make PDF file of graphs
plot.dmc(hss,hyper=TRUE,pll.chain=TRUE) 
plot.dmc(hss,hyper=TRUE,layout=c(2,4))
plot.dmc(hss,hyper=TRUE,layout=c(2,4),p.prior=pp.prior)
dev.off()
