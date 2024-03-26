# script to check multi-task data
# packages needed to be installed through terminal in the conda env 'R'
rm(list=ls())

library(dplyr)
library(ggplot2)
library(stringr)
library(gridGraphics)
library(grid)
library(gridExtra)
library(plotrix)
library(IDPmisc)

# define function to calculate standard error
SE <- function(x) {
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

# round dataframe columns function
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

# Set-up directories
if(Sys.getenv('HOME')=='/home/scotti') {
  ## means it is on a server
  main_dir <- "/home/scotti/2022/behaveModelling" # LISA
  data_dir <- file.path(main_dir, 'SST/rdatas')
  res_dir <- file.path(main_dir, 'SST/samples')
  pdf_dir <- file.path(main_dir, 'SST/pdfs')
  cores <- 16
} else if (Sys.getenv('HOME')=='/Users/scotti') { # running locally using mountainDuck
  main_dir <- '/Users/scotti/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/LISA/2022/behaveModelling'
  data_dir <- file.path(main_dir, 'SST/rdatas')
  res_dir <- file.path(main_dir, 'SST/samples')
  pdf_dir <- file.path(main_dir, 'SST/pdfs')  
  cores <- 4
}

#### SST
# General info:
# no. trials = 320 
# no. stop trials = 80
# 3s trials

setwd(main_dir)

sst_paths <- Sys.glob(file.path(main_dir, 'taskdata','*_task-SST_*_block-2_events.tsv'))

sst_dat <- c()
for (i in 1:length(sst_paths)){
  
  dat_tsv <- read.table(sst_paths[i], sep='\t', header=TRUE, stringsAsFactors = FALSE)
  if ('rt' %in% colnames(dat_tsv)){
    sst_dat <- rbind(sst_dat, dat_tsv)}
  
}

sst_dat <- sst_dat %>% select(trial_nr, response, event_type, phase, direction, stopsig_trial, # take only the interesting columns
                              subject, current_ssd, rt, choice_key, block_nr, null_trial)

sst_dat$subject <- as.factor(sst_dat$subject)
sst_dat$trial_nr <- sst_dat$trial_nr+1
sst_dat <- sst_dat[sst_dat$null_trial==0,]

#sst_dat$trial_nr <- sst_dat$trial_nr+(max(sst_dat$trial_nr,na.rm = TRUE)*(sst_dat$block_nr-1))

# extract one line per trial (this is probably overcomplicated... but it works.. )
sst_dat <- sst_dat[sst_dat$event_type == 'stimulus' | (sst_dat$event_type=='response' & !is.na(sst_dat$rt)),] # take only stimlus and reponse even types
ind=grep("stimulus",sst_dat$event_type) # index rows that are event type stimulus
sst_dat <- sst_dat[sort(c(ind[c(0,abs(diff(ind))==1)==1]-1,grep("response",sst_dat$event_type))),] # only keep stimulus rows when no response is given
sst_dat$response <- as.character(sst_dat$response)
sst_dat$response[sst_dat$event_type=='stimulus'] <- 'None'

sst_dat$type <- '' # Add trial types: succesful stop (SS), failed stop (FS), GO
sst_dat$type[sst_dat$stopsig_trial==1 & !sst_dat$choice_key %in% c('r','b')] <- 'SS'
sst_dat$type[sst_dat$stopsig_trial==1 & sst_dat$choice_key %in% c('r','b')] <- 'FS'
sst_dat$type[sst_dat$stopsig_trial==0] <- 'GO'

# add accuracy
sst_dat$correct <- 0
sst_dat$correct[sst_dat$choice_key=='r' & sst_dat$direction==1] <- 1
sst_dat$correct[sst_dat$choice_key=='b' & sst_dat$direction==0] <- 1

# change the type of some data
sst_dat$type <- as.factor(sst_dat$type)
sst_dat$direction <- as.factor(sst_dat$direction)
sst_dat$stopsig_trial <- as.factor(sst_dat$stopsig_trial)
sst_dat$current_ssd[sst_dat$current_ssd==-1] <- Inf

sst_dat <- droplevels(sst_dat)

sst_dat <- sst_dat %>% select(trial_nr, response, direction, stopsig_trial, subject, current_ssd, rt, type, correct, block_nr)

# Take a look at the structure of the data before we start analysis
str(sst_dat)
# 'data.frame':	1600 obs. of  9 variables:
#   $ trial_nr     : num  1 2 3 4 5 6 7 8 9 10 ...
# $ response     : Factor w/ 2 levels "m","z": 1 1 2 1 1 2 2 NA 1 1 ...
# $ direction    : Factor w/ 2 levels "0","1": 1 1 2 1 1 2 2 1 1 1 ...
# $ stopsig_trial: Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 2 1 1 ...
# $ subject      : Factor w/ 5 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ current_ssd  : num  Inf Inf Inf Inf Inf ...
# $ rt           : num  0.652 0.542 0.5 0.35 0.492 ...
# $ type         : Factor w/ 3 levels "FS","GO","SS": 2 2 2 2 2 2 2 3 2 2 ...
# $ correct      : num  1 1 1 1 1 1 1 0 1 1 ...

# create empty dataframe to append to
sst_res <- data.frame('sub'=0,'rt'=0,'acc'=0,'SS'=0,'SSperc'=0,'FS'=0,'FSperc'=0,'meanSSD'=0,'SSRT'=0,'gf'=0)

counter = 1
# Now we can loop through participants to calculate simple statistics
for (i in levels(sst_dat$subject)){
  
  this_sub <- sst_dat[sst_dat$subject == i,]    # data for this subject
  gdat <- this_sub[this_sub$stopsig_trial==0,]  # go trial data
  sdat <- this_sub[this_sub$stopsig_trial==1,]  # stop trial data
  mean_rt <- mean(gdat$rt,na.rm = TRUE)         # mean RT on go trials
  mean_acc <- mean(gdat$correct)                # mean accuracy on go trials
  num_ss <- nrow(sdat[sdat$type=='SS',])        # number of successful stops
  perc_ss <- (num_ss/nrow(sdat))*100            # percentage successful stops
  num_fs <- nrow(sdat[sdat$type=='FS',])        # number of failed stops
  perc_fs <- (num_fs/nrow(sdat))*100            # percentage of failed stops
  mean_ssd <- mean(NaRV.omit(this_sub$current_ssd))  # mean SSD
  ssrt <- median(this_sub$rt,na.rm = TRUE)-mean_ssd  # SSRT
  go_omis <- nrow(gdat[gdat$rt=='NA',])              # number of go omissions
  
  sst_res[counter,1:ncol(sst_res)] <- c(i,round(mean_rt,3),round(mean_acc,2),num_ss,perc_ss,
                                        num_fs,perc_fs,round(mean_ssd,2),round(ssrt,2),go_omis)
  
  counter=counter+1
}
sst_res

# START THE ANALYSIS

pdf(file.path(res_dir,'behavioural_summary_SST.pdf'),width = 13,height = 9) # make PDF file of graphs

group_dat <- data.frame('subject'=double(), 'totalTrials'=double(),'num_stop_trials'=double(),'num_go_trials'=double(),'stop_acc'=double(),'num_stops'=double(),'meanSSD'=double(),
                        'lastSSD'=double(), 'meanRT'=double(),'medianRT'=double(),'SSRT'=double())

# look at summary statistics per subject
for (i in levels(sst_dat$subject)){
  
  # separate go trial data from stop trial data
  all_dat <- sst_dat[sst_dat$subject==i,] # all data for this subject
  go_dat <- sst_dat[sst_dat$stopsig_trial==0 & sst_dat$subject==i,] # all go trial data for this subject
  stop_dat <- sst_dat[sst_dat$stopsig_trial==1 & sst_dat$subject==i,] # all stop trial data for this subject
  
  totalTrials <- nrow(all_dat) # total number of trials
  nstopTrials <- nrow(stop_dat) # number of stop trials
  pstopTrials <- (nstopTrials/totalTrials)*100 # percentage of trials that were stop trials
  ngoTrials <- nrow(go_dat) # number of go trials
  accStop <- nrow(stop_dat[stop_dat$response=='None',])/nrow(stop_dat) # % accuracy of stopping on stop trials
  numStop <- nrow(stop_dat[stop_dat$response=='None',]) # number of times the subject stopped on a stop trial
  meanSSD <- mean(stop_dat$current_ssd) # mean ssd 
  lastSSD <- tail(stop_dat$current_ssd, n=1) # last ssd 
  meanRT <- mean(go_dat$rt[go_dat$response!='None']) # mean reaction time 
  meanRT_left <- mean(go_dat$rt[go_dat$direction==0 & go_dat$response!='None']) # mean reaction time when stimulus was left
  meanRT_right <- mean(go_dat$rt[go_dat$direction==1 & go_dat$response!='None'])# mean reaction time when stimulus was right
  medianRT <- median(go_dat$rt[go_dat$response!='None'])#median reaction time
  accGo <- mean(go_dat$correct[go_dat$response!='None']) # accuracy on go trials when the subject responded
  accGoAll <- mean(go_dat$correct) # accuracy on go trials including trials where they didnt respond
  goOmissions <- length(go_dat$response[go_dat$response=='None']) # number of times they didnt respond on a go trial (should be low)
  SSRT <- medianRT - meanSSD # the stop signal reaction time
  
  # if too many stop or failed stops we will remove
  to_remove = 0
  if (accStop > 0.7 | accStop < 0.3){
    to_remove = 1
  }
  
  # print summary statistics for each subject
  cat("\n\n----------------------------\n\n", paste0('subject: ', i), paste0("\n\ntotal trials: ", totalTrials), paste0("\nnum stop trials: ", nstopTrials), paste0("\nnum go trials: ", ngoTrials), 
      paste0("\npercent stop trials (%): "), round(pstopTrials,2), paste0("\nstopping accuracy (%): ", round(accStop*100,2)), paste0("\nnumber of succesful stops: ", numStop), paste0("\nmean SSD (s): ", round(meanSSD,2)),  
      paste0("\nmean RT (s): ", round(meanRT,2)), paste0("\nlast SSD (s): ", round(lastSSD,2)),
      paste0("\nmean RT stimulus left (s): ", round(meanRT_left,2)), paste0("\nmean RT stimulus right (s): ", round(meanRT_right,2)), paste0("\nmedian RT (s): ", round(medianRT,2)), 
      paste0("\naccuracy go trials when responding (%): ", round(accGo*100,2)), paste0("\naccuracy all go trials (%): ", round(accGoAll*100,2)), paste0("\nGo omissions: ", goOmissions),
      paste0("\nSSRT (s): ", round(SSRT,2)), "\n\n")
  
  # append to dataframe
  group_dat <- rbind(group_dat, data.frame('subject'=i, 'totalTrials'=totalTrials, 'num_stop_trials'=nstopTrials, 'num_go_trials'=ngoTrials, 'stop_acc'=accStop, 'num_stops'=numStop, 
                                           'meanSSD'=meanSSD, 'lastSSD'=lastSSD, 'meanRT'=meanRT, 'medianRT'=medianRT, 'SSRT'=SSRT))
  
  # create some graphs
  par(mfrow=c(2,2))
  
  # plot density of RTs
  plot(density(go_dat$rt[go_dat$response!='None']), xlim=c(0,2), main=paste0("RT distribution subject ", i, ", SSRT: ", round(SSRT,2))
       , xlab='RT (s)')
  abline(v=c(meanSSD, meanRT, medianRT), lty=c(2,2,2), col=c('red','blue','green'))
  legend(1.5,max(density(go_dat$rt[go_dat$response!='None'])$x)-0.3,legend=c('ssd','mean RT','median RT'), col=c('red','blue','green'), lty=c(2,2,2))
  
  ##### NEW
  
  aplot = paste0('subject: ', i)
  bplot = paste0("total trials: ", totalTrials)
  cplot = paste0("num stop trials: ", nstopTrials)
  dplot = paste0("num go trials: ", ngoTrials)
  eplot = paste0("percent stop trials (%): ", round(pstopTrials,2))
  fplot = paste0("stopping accuracy (%): ", round(accStop*100,2))
  gplot = paste0("number of succesful stops: ", numStop)
  hplot = paste0("mean SSD (s): ", round(meanSSD,2))
  iplot = paste0("mean RT (s): ", round(meanRT,2))
  jplot = paste0("last SSD (s): ", round(lastSSD,2))
  kplot = paste0("mean RT stimulus left (s): ", round(meanRT_left,2))
  lplot = paste0("mean RT stimulus right (s): ", round(meanRT_right,2))
  mplot = paste0("median RT (s): ", round(medianRT,2))
  nplot = paste0("accuracy go trials w/ response (%): ", round(accGo*100,2))
  oplot = paste0("accuracy all go trials (%): ", round(accGoAll*100,2))
  pplot = paste0("Go omissions: ", goOmissions)
  qplot = paste0("SSRT (s): ", round(SSRT,2))
  
  plot(NA, xlim=c(0,11), ylim=c(0,11), bty='n',
       xaxt='n', yaxt='n', xlab='', ylab='')
  
  text(0,8,aplot, pos=4)
  text(0,7,bplot, pos=4)
  text(0,6,cplot, pos=4)
  text(0,5,dplot, pos=4)
  text(0,4,eplot, pos=4)
  text(0,3,fplot, pos=4)
  text(0,2,gplot, pos=4)
  text(0,1,hplot, pos=4)  
  
  text(5,8,iplot, pos=4)
  text(5,7,jplot, pos=4)
  text(5,6,kplot, pos=4)
  text(5,5,lplot, pos=4)
  text(5,4,mplot, pos=4)
  text(5,3,nplot, pos=4)
  text(5,2,oplot, pos=4)
  text(5,1,pplot, pos=4)
  
  #####
  
  # plot RT over course of the experiment
  plot(go_dat$trial_nr[go_dat$response!='None'], go_dat$rt[go_dat$response!='None'], 
       main=paste0("RT over experiment, subject ",i), xlab="Trial number", ylab="RT (s)", 
       ylim=c(0,max(sst_dat$rt[sst_dat$subject==i]+0.05,na.rm = TRUE)), type="l")
  
  # plot staircase
  plot(stop_dat$trial_nr, stop_dat$current_ssd, type="l", col="blue", ylim=c(0,max(NaRV.omit(sst_dat$current_ssd[sst_dat$subject==i]))+0.05),
       main=paste0("Staircase over experiment, subject ", i), xlab="Trial number", ylab="SSD (s)")
}

# round group table
group_dat[,2:11] <- round_df(group_dat[,2:11],2)
group_dat[nrow(group_dat)+1,2:11] <- round(colMeans(group_dat[,2:ncol(group_dat)]),2)

for (i in 1:ncol(group_dat[,2:11])){
  group_dat[nrow(group_dat),i+1] <- paste(round(mean(group_dat[,i+1]),2), '±', round(SE(group_dat[1:nrow(group_dat)-1,i+1]),2))
}

par(mfrow=c(1,1))

# plot &save summary of all participants with meaned column at bottom
plot(NA, xlim=c(0,8), ylim=c(0,8), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
addtable2plot(0 ,0,group_dat,bty="o",display.rownames=TRUE,hlines=TRUE,
              vlines=TRUE,title="Summary table\n")

dev.off()

# manually removing subs based on behavioural_summary_SST pdf
subs_to_remove = c(4,8, 19,27)

