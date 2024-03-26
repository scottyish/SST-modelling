# RUN ON PC

rm(list=ls())

# setwd('/home/scotti/Downloads/behaveModelling')
# res_dir <- '/home/scotti/Downloads/behaveModelling/SST2/rdatas'
# cores = 8

# Set-up directories
if(Sys.getenv('HOME')=='/home/scotti') {
  ## means it is on a server
  main_dir <- "/home/scotti/2022/behaveModelling" # LISA
  data_dir <- file.path(main_dir, 'SST/rdatas')
  res_dir <- file.path(main_dir, 'SST/samples')
  pdf_dir <- file.path(main_dir, 'SST/pdfs')
  cores <- 16
} else if (Sys.getenv('HOME')=='/Users/scotti') { # running locally using mountainDuck
  main_dir <-  '/Users/scotti/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/LISA/2022/behaveModelling'
  data_dir <- file.path(main_dir, 'SST/rdatas')
  res_dir <- file.path(main_dir, 'SST/samples')
  pdf_dir <- file.path(main_dir, 'SST/pdfs')  
  cores <- 4
}

setwd(main_dir)

source ("dmc/dmc.R")
load_model ("EXG-SS","exgSSprobit.R")
# tmp=load("SST2/rdatas/ss.RData")

load(file.path(res_dir,"ss_ntnu.RData"))

# STAGE 1
sss <- h.RUN.dmc(sss,cores=cores)
save(ss,sss,file=file.path(res_dir,"ss_stage1.RData"))

# Make hierarchical
p.prior <- sss[[1]]$p.prior
p1 <- get.p.vector(sss[[1]])[names(p.prior)]
s.prior <- prior.p.dmc(p1=p1,p2=p1,
                       dists=rep("gamma",length(p1)))
pp.prior=list(p.prior,s.prior)
hstart <- make.hstart(sss)
theta1 <- make.theta1(sss)
hss <- h.samples.dmc(nmc=100,p.prior,ss,pp.prior,
                     hstart.prior=hstart,theta1=theta1,thin=10)
rm(ss,sss); gc()

# STAGE 2
pp.prior <- attr(hss,"hyper")$pp.prior
hss  <- h.run.unstuck.dmc(hss, p.migrate = .05, cores = cores)
save(hss,file=file.path(res_dir,"ss_stage2.RData"))

hss <- h.run.converge.dmc(h.samples.dmc(nmc=100, samples=hss),
                          thorough=TRUE,nmc=50,cores=cores,finalrun=TRUE,finalI=250)
save(hss,file=file.path(res_dir,"ss_stage2.RData"))

hss <- h.run.dmc(h.samples.dmc(samples=hss,nmc=250),cores=cores,report=10)
save(hss,file=file.path(res_dir,"ss_stage2.RData"))

# POST PREDICT
ppss <- h.post.predict.dmc(hss,cores=cores)
save(hss,ppss,file=file.path(res_dir, "ss_post.RData"))

# SAVE CHAINS
pdf(file.path(pdf_dir,"ss_chains.pdf"),height=6,width = 8)
plot.dmc(hss,hyper=TRUE,pll.chain=TRUE) 
plot.dmc(hss,hyper=TRUE,layout=c(3,3))
plot.dmc(hss,hyper=TRUE,layout=c(3,3),p.prior=pp.prior)
dev.off()

# SAVE FITS
pdf(file.path(pdf_dir,"ss_fit.pdf"),height=6, width = 8)
plot.pp.dmc(ppss,model.legend = FALSE,layout=c(2,3))
dev.off()

# INFO CRITERIONS
h.IC.dmc(hss,DIC=TRUE)

gelman.diag.dmc(hss,hyper=TRUE)
gelman.diag.dmc(hss)

parsss <- summary.dmc(hss,hyper=TRUE)
round(parsss$quantiles[,c(1,3,5)],3)
parsss <- summary.dmc(hss)
save(hss,ppss,parsss,file=file.path(res_dir, "ss_post.RData"))

# SAVE IMAGE
save.image(file=file.path(res_dir,"img_save_ss.RData"))

# SAVE INDIVDUAL CHAINS
dir.create(file.path(pdf_dir, paste0("chains_pp")))
for (i in 1:length(hss)){
  pdf(file.path(pdf_dir, paste0("chains_pp/ss_chains_sub_",names(hss)[i],".pdf")),height=6,width=8)
  plot.dmc(hss[i],hyper=FALSE,pll.chain=TRUE)
  plot.dmc(hss[i],hyper=FALSE,layout=c(3,3))
  plot.dmc(hss[i],hyper=FALSE,layout=c(3,3),p.prior=p.prior)
  dev.off()
}

# SAVE INDIVIDUAL FITS
pdf(file.path(pdf_dir,"ss_fit_pp.pdf"),height=12,width = 16)
for (i in 1:length(ppss)){
  plot.pp.dmc(ppss[[i]], style='cdf',layout=c(2,2),ylim=c(0,1),fits.pcol="grey",
              model.legend=TRUE,dname=paste("Sub",names(hss)[i],"data"),mname=paste("Sub",names(hss)[i],"model"), aname="")
  
}
dev.off()

# SAVE IMAGE LAST TIME
save.image(file=file.path(res_dir,"img_save_ss.RData"))
save.image(file="img_save_ss.RData")