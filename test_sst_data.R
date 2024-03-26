# check and sort/analyze SST data
# This script find outliers and checks the data
# the output from this script should be passed to ss_mod_setup.R

rm(list=ls())

library(gamlss.dist)

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

setwd(main_dir)

# setwd('/Users/scotti/surfdrive/Projects/7T_SST_MSIT_NTNU')

source("dmc/dmc.R")
load_model ("EXG-SS","exgSS.R")

load(file.path(data_dir,'dat_for_test_data.RData'))

########## DATA STRUCTURE #########

head(dat)
### FINAL
# R     S s  SSD        RT SS
# 11 RIGHT right 2 0.20 0.6176246 SS
# 21 RIGHT right 2  Inf 0.5684693 GO
# 30  LEFT  left 2  Inf 0.3685658 GO
# 41 RIGHT right 2  Inf 0.4989910 GO
# 50    NR right 2 0.15        NA SS
# 60  LEFT  left 2 0.20 0.5483389 SS

####### TEST ALL SUBJECTS FOR EXCLUSIONS AND OUTLIERS ######
dat$R <- factor(dat$R, levels = c("NR","LEFT","RIGHT")) # put factor levels in correct order
#sort_labs <- paste(sort(as.integer(levels(dat$s))))
#dat$s <- factor(dat$s, levels=sort_labs) # fix lexicographic sorting of numbers

remove_subs <- c() # create empty vector of subs to remove

str(dat)
## FINAL
# 'data.frame':	7100 obs. of  6 variables:
#   $ R  : Factor w/ 3 levels "NR","LEFT","RIGHT": 3 3 2 3 1 2 2 1 1 3 ...
# $ S  : Factor w/ 2 levels "left","right": 2 2 1 2 2 1 1 2 1 2 ...
# $ s  : Factor w/ 37 levels "2","3","4","5",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ SSD: num  0.2 Inf Inf Inf 0.15 ...
# $ RT : num  0.618 0.568 0.369 0.499 NA ...
# $ SS : Factor w/ 2 levels "GO","SS": 2 1 1 1 2 2 1 2 2 1 ...

lapply(dat,levels)
### NEW
# $R
# [1] "NR"    "LEFT"  "RIGHT"
# 
# $S
# [1] "left"  "right"
# 
# $s
# [1] "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "23" "24"
# [23] "25" "26" "27" "29" "31" "32" "33" "34" "35" "37" "38" "39" "40" "41" "42"
# 
# $SSD
# NULL
# 
# $RT
# NULL
# 
# $SS
# [1] "GO" "SS"

# starting with 37 subjects

dat <- dat[dat$RT >= 0.15 | is.na(dat$RT)==1,]   # Remove RTs < 200ms as these could not have be conscious decisions
dat <- dat[dat$RT <= 1.6 | is.na(dat$RT)==1,]
is.stop <- is.finite(dat$SSD) # find stop signal trials
gdat <- dat[!is.stop,]        # find go trial data
ssdat <- dat[is.stop,]        # find stop trial data

# Go omissions
is.nr <- gdat$R=="NR"         # which trials were go omissions on go trials 
# get rid of people with > 10 go omissions
sort(round(tapply(is.nr,gdat$s,mean)*100,2)) # percentage go omissions
# 2     6    20    23    31    32    34    42     3     5    25    33    39    16    18    21    24    37    40 
# 0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.67  0.67  0.67  0.67  0.67  1.33  1.33  1.33  1.33  1.33  1.33 
# 41     7    12     9    11    15    17    26    29    14    35    38     4    10     8    19    13    27 
# 1.33  2.00  2.00  2.67  2.67  3.33  4.67  4.67  5.33  6.67  6.67  6.67  8.00  8.67 10.00 15.33 22.67 28.67 

sort(tapply(is.nr,gdat$s,sum)) # sum of go omissions
# 2  6 20 23 31 32 34 42  3  5 21 25 33 39 40 16 18 24 37 41  7 12  9 11 15 17 26 29 14 35 38  4 10  8 13 19 27 
# 0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  2  2  2  2  2  3  3  4  4  5  7  7  8 10 10 10 12 13 15 17 23 43 

remove_subs <- append(remove_subs, c(4,10,8,19,13,27)) # remove subs with > 10 go omissions

# Go accuracy on trials where subjects responded
crct <- gdat$R[gdat$R!='NR'] == toupper(gdat$S[gdat$R!='NR'])
round(sort(tapply(crct,gdat[gdat$R!='NR',"s"],mean)),3)
# 18    23     8    29    11    16     5    33    39     2     3     4     6     7     9    10    12    13    14 
# 0.980 0.980 0.993 0.993 0.993 0.993 0.993 0.993 0.993 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 
# 15    17    19    20    21    24    25    26    27    31    32    34    35    37    38    40    41    42 
# 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000

remove_subs <- append(remove_subs, c()) # remove subs with < 95% accuracy

# check race model 
sort(round(tapply(gdat$RT,gdat$s,mean,na.rm=TRUE)-tapply(ssdat$RT,ssdat$s,mean,na.rm=TRUE),3))
# 13     27     23     34     16     10     19     32      8     39     35      6     20     33      7     11 
# -0.073 -0.024  0.011  0.019  0.028  0.029  0.032  0.047  0.052  0.054  0.056  0.059  0.071  0.075  0.076  0.079 
# 2     29     12     21     42      5     24     17     41     31      9     18     40      3     14     38 
# 0.082  0.083  0.084  0.086  0.087  0.092  0.095  0.097  0.098  0.103  0.104  0.104  0.117  0.119  0.120  0.121 
# 37     25     26     15      4 
# 0.125  0.137  0.139  0.140  0.176 

remove_subs <- append(remove_subs, c(13, 27)) # remove subs with results < 0 (means that ss RTs were longer than go RTs)

# check stopping accuracy
resp <- ssdat$R =='NR'
round(sort(tapply(resp, ssdat[,"s"], mean)),2)
# 23   27   13   32   39   16   33   34   41    6    9   18   40   42    2    3   11   29   31    5   20   21   25 
# 0.42 0.44 0.48 0.48 0.48 0.50 0.50 0.50 0.50 0.52 0.52 0.52 0.52 0.52 0.54 0.54 0.54 0.54 0.54 0.56 0.56 0.56 0.56 
# 35   37   38   15    4   10   12   17   24    7    8   14   26   19 
# 0.56 0.56 0.56 0.58 0.60 0.60 0.60 0.60 0.60 0.62 0.62 0.62 0.62 0.66 

remove_subs <- append(remove_subs, c(19)) # remove subs either < 35% or > 65% accurate

###### EXTRA CHECKS #######

# plot all ihibition functions per participant
par(mfrow=c(2,3))
for (i in levels(dat$s)){
  tapply(!is.na(dat$RT[dat$s==i]),dat[dat$s==i,c("SS","SSD")],mean)
  plot_SS_if.dmc(dat[dat$s==i,], main = paste0("Inhibition function, sub - ", i))
}

# plot all staircases per participant
for (i in levels(dat$s)){
  plot(1:length(dat$SSD[dat$SS=="SS"&dat$s==i]), dat$SSD[dat$SS=="SS"&dat$s==i], type="l", main=paste0("sub - ", i, ", staircase over exp"))
}

# RTs over experiment per participant
for (i in levels(dat$s)){
  plot(1:length(dat$RT[dat$SS=="GO"&dat$s==i]), dat$RT[dat$SS=="GO"&dat$s==i], type="l", main=paste0("sub - ", i, ", RT over exp"))
}

# remove_subs <- append(remove_subs, c(40,43,52,56, 60, 13,22,38))

#####################################################################
######### REMOVE EXCLUSIONS FROM THE REST OF THE ANALYSIS ###########
dat <- dat[!dat$s %in% remove_subs,]
dat <- droplevels(dat)
is.stop <- is.finite(dat$SSD) # find stop signal trials
gdat <- dat[!is.stop,]        # find go trial data
ssdat <- dat[is.stop,]        # find stop trial data
#####################################################################
#####################################################################

par(mfrow=c(1,1))

# Go RT
head(sort(gdat$RT))
tail(sort(gdat$RT))
hist(gdat$RT,breaks="fd")

# Stop RT
head(sort(ssdat$RT))
tail(sort(ssdat$RT))
hist(ssdat$RT,breaks="fd")

# Show the different SSDs:
sort(tapply(as.character(dat$SSD),dat[,c("SS")],unique)$SS)
# [1] "0.05" "0.1"  "0.15" "0.2"  "0.25" "0.3"  "0.35" "0.4"  "0.45" "0.5"  "0.55" "0.6"  "0.65" "0.7"  "0.75" "0.8" 

# Show the number of trials for each SSD:
Ns = tapply(dat$RT,dat$SSD,length)
Ns
# 0.05  0.1 0.15  0.2 0.25  0.3 0.35  0.4 0.45  0.5 0.55  0.6 0.65  0.7 0.75  0.8  Inf 
# 27   64  136  201  178  147  130  137  132   97   75   57   48   44   24    3 4499 

# Show response rate:
tapply(!is.na(dat$RT),dat[,c("SS")],mean)
# GO        SS 
# 0.9804401 0.4566667 

# Response rate broken down by SSD & corresponding inhibition function:
tapply(!is.na(dat$RT),dat[,c("SS","SSD")],mean)
plot_SS_if.dmc(dat)  #P(Respond) increases as a function of SSD, as it should

# Plot median signal-respond RT per SSD:
tapply(dat$RT,dat$SSD,median,na.rm=TRUE) 
plot_SS_srrt.dmc(dat) # Median SRRT increases as a function of SSD, as it should

# Show number of signal-respond RTs per SSD:
Nr = tapply(!is.na(dat$RT),dat[,c("SS","SSD")],sum)[2,]
Nr
# 0.05  0.1 0.15  0.2 0.25  0.3 0.35  0.4 0.45  0.5 0.55  0.6 0.65  0.7 0.75  0.8  Inf 
# 6   21   46   93   83   74   55   58   66   53   37   31   19   22   18    3   NA 

# Signal-respond RTs (stop dist, red line) should be faster than go RTs:
hist(dat$RT[dat$SS=="GO"],breaks="fd",main="SRRT vs. Go RT",freq=F,ylim=c(0,8)) # Go dist
lines(density(dat$RT[dat$SS=="SS"],na.rm=T),col="red",lwd=2) # Stop dist

# test the racce model
test_race <- function(s,ssdat,gdat) {
  grt <- gdat[gdat$s==s,"RT"]
  srt <- ssdat[ssdat$s==s,"RT"]
  t.res <- t.test(grt[!is.na(grt)],srt[!is.na(srt)])
  if (t.res$p.value >= 0.05){
  print(i)
  print(t.res$p.value)
  }
  # print(t.test(grt[!is.na(grt)],srt[!is.na(srt)]))
}

for (i in levels(gdat$s)){
test_race(i,ssdat,gdat)
}
  
unique(remove_subs)
# [1]  4 10  8 19 13 27

save(dat, file=file.path(data_dir, 'prepped_data_sst.RData'))


#   dat$correct <- as.numeric(dat$S)==as.numeric(dat$R)-1 # make correct oclumn for ease



