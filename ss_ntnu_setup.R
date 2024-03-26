# NEW MOD SETUP
rm(list=ls())

# run locally to setup model and then ss_ntnu_run.R on LISA

# Set-up directories
if(Sys.getenv('HOME')=='/home/scotti') {
  ## means it is on a server
  main_dir <- "/home/scotti/2022/behaveModelling" # LISA
  data_dir <- file.path(main_dir, 'SST/rdatas')
  res_dir <- file.path(main_dir, 'SST/samples')
  pdf_dir <- file.path(main_dir, 'SST/pdfs')
  cores <- 16
} else if (Sys.getenv('HOME')=='/Users/scotti') { # running locally using mountainDuck
  main_dir <- '/Users/scotti/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/LISA/home/scotti/2022/behaveModelling'
  data_dir <- file.path(main_dir, 'SST/rdatas')
  res_dir <- file.path(main_dir, 'SST/samples')
  pdf_dir <- file.path(main_dir, 'SST/pdfs')  
  cores <- 4
}

setwd(main_dir)

# libraries
library(dplyr)

# load dat shit
load(file.path(data_dir, 'prepped_data_sst.RData'))

# Fitting
source ("dmc/dmc.R")
load_model ("EXG-SS","exgSSprobit.R")

model <- model.dmc(
  # SS stands for trial type (GO or Stop-signal [SS]):
  factors=list(S=c("left","right"),SS=c("GO","SS")),
  # NR stands for "No response", i.e., go omission & successful inhibitions:
  responses=c("NR","LEFT","RIGHT"),
  # Match scores correct responses for each GO stimulus as usual, and also scores
  # the "correct" stimulus corresponding to an NR response, but the latter has
  # no effect (it just avoids a standard check making sure that each response is scored):
  match.map=list(M=list(left="LEFT",right="RIGHT",left="NR")),
  p.map=list(mu="M",sigma="M",tau="M",muS="1",sigmaS="1",tauS="1",tf="1",gf="1"),
  # no go failures (low go omissions in data)
  # constants=c(gf=0.0001), # if probit scale this causes 50% gfs
  constants = c(mu.false=1e6,sigma.false=.001,tau.false=.001),
  type="exgss")

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

# plor & save priors
pdf(file.path(pdf_dir,"ss_priors.pdf"),height=6,width = 8)
par(mfcol=c(2,6)); for (i in names(p.prior)) plot.prior(i,p.prior)
dev.off()

# check some things before running
length(attr(model,"p.vector"))
head(dat,50)
str(dat)

### The real thing
# Generate start points
ss <- data.model.dmc(dat,model)
sss <- h.samples.dmc(nmc=100,p.prior,ss)
save(ss,sss,file=file.path(res_dir,"ss_ntnu.RData"))






## OLD


# # start sampling
# samples <- samples.dmc(nmc=100,p.prior,data)
# samples <- run.dmc(samples,report=1,cores=cores,p.migrate=.05)
# plot.dmc(samples,pll.chain=TRUE)
# 
# samples1 <- RUN.dmc(samples,cores=cores,report=1,max.try=100,verbose=TRUE)
# 
# plot.dmc(samples1,layout=c(2,6))
# plot.dmc(samples1,pll.chain=TRUE)
# 
# gelman.diag.dmc(samples1)
# 
# tmp=check.recovery.dmc(samples1,p.vector)
# # Hard to judge tf and gf on the probit scale, so put on probability
# tmp[,c("tf","gf")] <- pnorm(tmp[,c("tf","gf")])
# tmp["Median-True",] <- tmp["50% Estimate",]-tmp["True",] 
# attributes(tmp)$ci50 <- NULL
# round(tmp,3)



# tmp <- load("ss.RData")



