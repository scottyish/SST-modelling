# create RData file for test_sst_dat.R

rm(list=ls())

library(dplyr)
library(ggplot2)
library(stringr)
library(gridGraphics)
library(grid)
library(gridExtra)
library(plotrix)
library(IDPmisc)

options(max.print=1000000)

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

# main_dir <- '/Users/scotti/surfdrive/Projects/7T_SST_MSIT_NTNU' # where is the data?
# res_dir <- '/Users/scotti/surfdrive/Projects/7T_SST_MSIT_NTNU/data/sstdata'

# LOAD DATA
# sst_dat <- read.csv(file.path(mainDir, 'DataSST.csv'), sep='\t')

sst_paths <- Sys.glob(file.path(main_dir, 'SST','taskdata','*_task-SST_*_block-2_events.tsv'))

sst_dat <- c()
for (i in 1:length(sst_paths)){
  
  dat_tsv <- read.table(sst_paths[i], sep='\t', header=TRUE, stringsAsFactors = FALSE)
  if ('rt' %in% colnames(dat_tsv)){
    sst_dat <- rbind(sst_dat, dat_tsv)}
  
}

# take only the interesting columns
# sst_dat <- sst_dat %>% dplyr::select(trial_nr, key_press, direction, stop_trial, correct, timestamp,
#                               subjectNumber, st_sig_delay, rt, block, response_button_left, response_button_right)

sst_dat <- sst_dat %>% dplyr::select(trial_nr, response, event_type, phase, direction, stopsig_trial, # take only the interesting columns
                              subject, current_ssd, rt, choice_key, block_nr, null_trial)

# extract one line per trial (this is probably overcomplicated... but it works.. )
sst_dat <- sst_dat[sst_dat$event_type == 'stimulus' | (sst_dat$event_type=='response' & !is.na(sst_dat$rt)),] # take only stimlus and reponse even types
ind=grep("stimulus",sst_dat$event_type) # index rows that are event type stimulus
sst_dat <- sst_dat[sort(c(ind[c(0,abs(diff(ind))==1)==1]-1,grep("response",sst_dat$event_type))),] # only keep stimulus rows when no response is given
sst_dat$response <- as.character(sst_dat$response)
sst_dat$response[sst_dat$event_type=='stimulus'] <- 'None'

#sort some shit
sst_dat$subject <- as.factor(sst_dat$subject)
sst_dat$trial_nr <- sst_dat$trial_nr+1
sst_dat <- sst_dat[sst_dat$null_trial==0,]
sst_dat$choice_key <- as.character(sst_dat$choice_key)

# change response keys
sst_dat$choice_key[sst_dat$choice_key == 'b'] <-'LEFT'
sst_dat$choice_key[sst_dat$choice_key == 'r'] <-'RIGHT'
sst_dat$choice_key[!sst_dat$choice_key %in% c('LEFT','RIGHT')] <- 'NR'

# Change direction of arrow text
sst_dat$direction[sst_dat$direction==0] <- 'left'
sst_dat$direction[sst_dat$direction==1] <- 'right'

# Add tiral types: succesful stop (SS), failed stop (FS), GO
sst_dat$type <- ''
sst_dat$type[sst_dat$stopsig_trial==1 & !sst_dat$choice_key %in% c('LEFT','RIGHT')] <- 'SS' # sucessful stop
sst_dat$type[sst_dat$stopsig_trial==1 & sst_dat$choice_key %in% c('LEFT','RIGHT')] <- 'SS' # failed stop, but SST modellnig just needs SS to denote stop trial
sst_dat$type[sst_dat$stopsig_trial==0] <- 'GO'

# change the type of some data
sst_dat$type <- as.factor(sst_dat$type)
sst_dat$direction <- as.factor(sst_dat$direction)
sst_dat$stopsig_trial <- as.factor(sst_dat$stopsig_trial)
sst_dat$current_ssd[sst_dat$stopsig_trial==0] <- Inf

# Remove stuff
sst_dat <- droplevels(sst_dat)

# Take only interesting columns
#sst_dat <- sst_dat %>% select(trial_nr, key_press, direction, stop_trial, subjectNumber, st_sig_delay, rt, type, correct, session)
sst_dat <- sst_dat %>% dplyr::select(choice_key, direction, subject, current_ssd, rt, type)

sst_dat <- sst_dat %>% rename(S = direction, SS = type, s = subject, RT = rt, R = choice_key, SSD = current_ssd)

# last things
sst_dat$R <- as.factor(sst_dat$R)

# Take a look at the structure of the data before we start analysis
str(sst_dat)
dat <- sst_dat

save(dat, file= file.path(res_dir, 'dat_for_test_data.RData'))