#### FLOCKTURE vs STRUCTURE


#### load some useful libraries ####
library(ggplot2)
library(stringr)
library(plyr)
library(reshape2)
library(miscFuncs)
library(dplyr)  # for the %>% operator
library(xtable)

##### Top-level settings and paths and stuff  ####

### Set paths
# Your R working directory should be the flock-compare directory
# and flockture and slg_pipe should be directories that sit beside flock-compare
# one directory level above.
# If you are using RStudio and opened the flock-compare project you should
# be good to go.

# now, we want to get absolute paths to a lot of the directories we will be using:
flockcommentDIR <- normalizePath("../flock-compare") # main working directory
flocktureDIR <- normalizePath("../flockture") #  this will be where we simulate all of our datasets
slg_pipeDIR <- normalizePath("../slg_pipe") # this is where clump and distruct live
StructBinaryPath <- normalizePath("../slg_pipe/inputs/ForStructure/structure2.3.4")

# now, this is a bit of a hack.  We need some of the binaries to be in a bin directory
# in the clump_and_distruct directory.  So we will just copy everything in the StructureArea bin to there.
file.copy(file.path(slg_pipeDIR, "inputs/ForStructure/bin"),
          file.path(slg_pipeDIR, "inputs/ForStructure/clump_and_distruct/"),
          recursive = TRUE
)
# and, another issue, the distruct executable in there does not see to run
# on Yosemite.  We need distruct_BIG anyway, probably, so copy that over.
file.copy(file.path(slg_pipeDIR, "inputs/ForStructure/clump_and_distruct/bin/distruct_BIG"),
          file.path(slg_pipeDIR, "inputs/ForStructure/clump_and_distruct/bin/distruct"),
          overwrite = TRUE
)
# Also, we have have to make this directory to put N.txt into it later during clumping and distructing.
dir.create(file.path(slg_pipeDIR, "/inputs/ForStructure/input"))  
 


### Set genetic marker type to simulate
marker=2 # RUN WITH MICROSATELLITES (=1) OR SNPS (=2) 

### Other options to set
qi_indLoss <- F # set true if you want graphs of qi values and loss by simulation

##### PART 1: Simulate data sets ####

### Set up the parameters controlling the simulations 
Seedset <-  3
Reps <- 9 # how many reps do we want to run for each program
Npops <- 5
N <- 2000

## set up migration rates to use 
if(marker == 1) {
  MGrate <- c(20,24,28,30,31,32,36,37,37.5,38,40,50,60,75,90,125)
} else if(marker == 2) {
  MGrate<-c(20,24,28,30,31,34,36,36.8,37.5,38,40,50,60,75,80,87)
} else {
  stop("Unrecognized value for variable marker")
}

DatNum <- length(MGrate) # how many datasets do we have?


### Simulate some datasets and store the input files in Dat* folder

# create directories to for simulation output
setwd(flockcommentDIR)
if(any(lapply(paste("SimDat", 1:(DatNum*Seedset), sep = ""), dir.create) == FALSE)) {
  stop("Failed to create SimDat directory. Perhaps they already exist? If so, put 'em somewhere else!")
}

if(marker==1) {

  ## FOR uSATs
  # loop over our MGrate vector to create simulated datasets and store them in 
  # folders named SimDatX where X is the number of the simulation
  # details on the simulation are stored in SimDets to see the pairwise Fst Values etc.
  Nloci<-15
  
  setwd(file.path(flockcommentDIR, "simdata"))
  
  seedsms<-readLines("uSat_seedsms.txt")
  seedsms2geno<-readLines("uSat_seedsms2geno.txt")
  
  i <- 0  # to number data sets
  for (s in 1:Seedset) {
    # set seeds
    cat(seedsms[s], file = "seedms")
    cat(seedsms2geno[s], file = "ms2geno_seeds")
     
    for(x in 1:DatNum) { 
      
      # run the shell command to generate baseline and structure file
      message(paste("Simulating microsat data set x =", x, "and s=", s))
      system(paste("./sim_data_MIG_RATE.sh", MGrate[x], "> SimDets.txt"))
      
      # compute average pairwise Fst between all the simulated pops
      tmp <- readLines("SimDets.txt")  # read in the simulation log file
      tmp <- tmp[str_detect(tmp, "^PAIRWISE_FST")]  # get the lines that have pairwise Fst on them
      strsplit(tmp, ":") %>%             # split on the colon
        sapply(., "[", 4)  %>%               # grab the fourth field
        as.numeric  %>%                      # make it numeric
        mean  %>%                            # compute mean
        cat(., file = "AvgFst.txt", sep = "\n")          # write it to file
      
       
      # add pop label to the structure input file and save it into file SimDatInX
      # here is some fun regex foo to do that  
      i <- i + 1
      tmp <- readLines("struct_input_1.txt")
      str_match(tmp, "^(Pop_([0-9]+)_BaseInd_[0-9]+)(.*$)")[, -1] %>%    # pick out and put back that number
        write.table(., sep = "  ", quote = FALSE, row.names = FALSE, col.names = FALSE, 
                    file = paste("SimDatIn", i, sep = ""))
      
        
      # now, name some files we wish to move into "../SimDatX"
      files2mov<-c("AvgFst.txt", 
                   "SimDets.txt", 
                   dir(pattern = "BaseFile*"),    #use * just in case we want to generate more inputfiles
                   dir(pattern = "struct_input*"), 
                   paste("SimDatIn", i, sep = "")
                   ) 
      
      lapply(files2mov, function(x) file.rename(x, paste("../SimDat", i, "/", x, sep = "")))
    }
  }#do for uSats
}

if(marker==2) {
  Nloci<-96
  setwd(paste(flockcommentDIR,"/simdata",sep=""))
  seedsms<-readLines('SNPs_seedsms.txt',warn=F)
  seedsms2geno<-readLines('SNPs_seedsms2geno.txt',warn=F)
  index<-seq(from=1,to=DatNum*Seedset,by=DatNum)
  
  i <- 0  # to number data sets
  for (s in 1:Seedset) {
    # set seeds
    cat(seedsms[s], file = "seedms")
    cat(seedsms2geno[s], file = "ms2geno_seeds")
    
    for(x in 1:DatNum){ 
      
      # run the shell command to generate baseline and structure file
      message(paste("Simulating SNP data set x =", x, "and s=", s))
      system(paste("./sim_data_snps_MIG_RATE.sh", MGrate[x], "> SimDets.txt"))
      
      #what is the avg. pariwise FST?
      # compute average pairwise Fst between all the simulated pops
      tmp <- readLines("SimDets.txt")  # read in the simulation log file
      tmp <- tmp[str_detect(tmp, "^PAIRWISE_FST")]  # get the lines that have pairwise Fst on them
      strsplit(tmp, ":") %>%             # split on the colon
        sapply(., "[", 4)  %>%               # grab the fourth field
        as.numeric  %>%                      # make it numeric
        mean  %>%                            # compute mean
        cat(., file = "AvgFst.txt", sep = "\n")          # write it to file
        
      
      # add pop label to the structure input file and save it into file SimDatInX
      # here is some fun regex foo to do that  
      i <- i + 1
      tmp <- readLines("struct_input_1.txt")
      str_match(tmp, "^(Pop_([0-9]+)_BaseInd_[0-9]+)(.*$)")[, -1] %>%    # pick out and put back that number
        write.table(., sep = "  ", quote = FALSE, row.names = FALSE, col.names = FALSE, 
                    file = paste("SimDatIn", i, sep = ""))
      
        
      # now, name some files we wish to move into "../SimDatX"
      files2mov<-c("AvgFst.txt", 
                   "SimDets.txt", 
                   dir(pattern = "BaseFile*"),    #use * just in case we want to generate more inputfiles
                   dir(pattern = "struct_input*"), 
                   paste("SimDatIn", i, sep = "")
                   ) 
      
      lapply(files2mov, function(x) file.rename(x, paste("../SimDat", i, "/", x, sep = "")))
    }
  }
}  # end if(marker == 2) 

##### PART 2: Run flockture on all of these datasets #####
setwd(flocktureDIR)

# source the flockture code... it should be compiled first 
# by running 'make' in terminal within scr folder
source("R/flockture.R")

### Run FLOCKTURE 
# We ran FLOCK 9 times with 50 reps and 20 interations...
datasets<-1 #how many datasets for each of the simulations did we make?
for (Ds in 1:(DatNum*Seedset)){
  
  dat<-paste(flockcommentDIR,'/SimDat',Ds,'/struct_input_',datasets,'.txt',sep='')
  
  # read in a data set:
  D <- read.table(dat, row.names = 1)
  
  # define some starting conditions if you want:
  #sc <- c(c(rep(2,times=14),rep(1,times=15),rep(2,times=10),rep(1,times=11))
  #)
  
  FLCTR_PlateauRec<-vector()
  Flockture.Runtime<-vector()
  FlocktureReps<-Reps #how many times are we running flockture
  for (P in 1:FlocktureReps){
    message("Doing Flockture run on SimDat ", Ds, " out of ", DatNum * Seedset, "    Rep ", P, " of ", FlocktureReps)
    # run flockture on it and grab the results out
    ptm<-proc.time()
    catch <- run_flockture_bin(D, K = 5, iter = 20, reps = 50)
    out <- slurp_flockture_dumpola()
    
    # summarize by plateaus and add that to our output
    psum <- plateau_summarize(out)
    
    # now, draw lines of log_prob against iterations for the 
    # reps, and color them according to plateau length
    ggplot(psum$log_probs, aes(x = iter, y = log.prob, group = rep, color = factor(plat.len))) + 
      geom_line() +
      scale_color_discrete(name="Plateau Length")
    Flockture.Runtime[P]<-(proc.time()-ptm)[3]
    
    # here is a vector of the different plateau lengths in sorted order:
    sapply(psum$plateau_list, length)
    FLCTR_PlateauRec1<-vector()
    FLCTR_PlateauRec1<-sapply(psum$plateau_list, length)
    FLCTR_PlateauRec[P]<-paste(FLCTR_PlateauRec1,collapse=",")
    
    # Which FLOCKTURE rep is the best?
    # We should take the result with the longest plateau, 
    # but if all plateaus are the same length we should take the one wiht the highest log.prob
    # Longest plateau does not mean highest log.prob! So we need to do one and then the next.
    
    #which reps belong to the largest plateau?
    Max.pl<-max(psum$log_probs[psum$log_probs$iter==20,5])
    Max.reps<-psum$log_probs[psum$log_probs$iter==20,][psum$log_probs[psum$log_probs$iter==20,5]==Max.pl,]
    Best.rep<-Max.reps[which(Max.reps$log.prob==max(Max.reps$log.prob))[1],1]
    
    Fkture.ref<-out$indivs[out$indivs$rep==Best.rep,4]
    
    #### Write out qi values 
    ResOut<-round(out$indivs[out$indivs$rep==Best.rep,(5:9)],3)    
    write.table(format(ResOut,digits=4),paste('FLOCKTUREqi_REP',P,'.txt',sep=""),sep=' ',row.names=FALSE,col.names=FALSE,quote=FALSE)
    
    #paste results into StructureOutput
    #extract results
    #pasteCmd<-paste("awk 'NR==FNR{a[FNR]=$0;next} {print a[FNR],$0}' ",flockcommentDIR,"/Flockture2ClmpDstrct_Skeleton  ~/flockture/",paste('FLOCKTUREqi_REP',F,'.txt',sep="")," > ", flockcommentDIR,"/SimDat",Ds,'/StructOuput_genos_slg_pipe.txt_dat00',Ds,'_k005_Rep03',F,'.txt_f', sep="")
    pasteCmd<-paste("awk 'NR==FNR{a[FNR]=$0;next} {print a[FNR],$0}' ",flockcommentDIR,"/Flockture2ClmpDstrct_Skeleton  ",flocktureDIR,paste('/FLOCKTUREqi_REP',P,'.txt',sep="")," > ", flockcommentDIR,"/SimDatIntermediate", sep="")
    system(pasteCmd)
    
    #paste into the STRUCTURE output skeleton
    File<-paste(flockcommentDIR,"/SimDatIntermediate",sep="") 
    
    sedcmd<- paste("sed '/INSERT_INDQ_HERE/ r ",File, "' <",flockcommentDIR,"/StructClmpDstrct_skeleton >",flockcommentDIR,"/Intrmd.txt",sep="")
    system(sedcmd)
    
    sedcmd<-paste("sed '/INSERT_INDQ_HERE/d' '",flockcommentDIR,"/Intrmd.txt' >",flockcommentDIR,"/SimDat",Ds,"/StructOuput_genos_slg_pipe.txt_dat00",Ds,"_k005_Rep03",P,".txt_f",sep="")
    system(sedcmd)
    
    #### MOVE or DELETE the original qi files from flockture
    # move the flockture results into the SimDat folder
    cmd<-paste('mv FLOCKTUREqi* ', flockcommentDIR,'/SimDat',Ds,sep="" )
    system(cmd)
    
  }#over P reps of Flockture
  
  write.table(x=FLCTR_PlateauRec,file='FlockturePlateaus.csv',sep=",",quote=FALSE,append=TRUE,row.names=FALSE,col.names=FALSE)
  write.table(x=Flockture.Runtime,file='FlocktureRuntime.csv',sep=",",quote=FALSE,append=TRUE,row.names=FALSE,col.names=FALSE)
  
  # move all the stuff 
  files2mov<-c('Flockture*.csv')
  length(files2mov)
  path<-lapply(1:length(files2mov),function (y) paste('mv ',files2mov[y],' ',flockcommentDIR,'/SimDat',Ds,sep=""))
  lapply(1:length(files2mov),function(y) system(path[[y]][1]))
}#over Ds datasets

##### PART 3:  Now run STRUCTURE ########################
if(marker==1){
  mark<-'uSat'
}
if(marker==2){
  mark<-'SNPs'
}


main_param<-paste(flockcommentDIR,'/mainparams_',sep="") #note this has an _ because it will change depending on the model that is run
ex_param<-paste(flockcommentDIR,'/extraparams',sep="")#for these models we don't have any extraparams so we can use the same blank file

# need to run 3 models:  AdmixCorr, NonAdmixCorr, NonAdmixNonCorr - might want to make this user specified
# seach SimDat file Rep times
for (SDt in 1:(DatNum*Seedset)){ # goes into different SimDat folders
  # lets make some vectors to store some processing times
  AdmixCorrRuntime<-vector()
  NoAdmixCorrRuntime<-vector()
  NoAdmixNonCorrRuntime<-vector()
  for (R in 1:Reps){
    message(paste("Running Structure:  SimDat", SDt, "of", DatNum*Seedset, "   Rep Number", R, "of", Reps))
    
    #AdmixCorr
    message("     AdmixCorr")
    ptm<-proc.time()
    system(paste(StructBinaryPath,' -m ',main_param,'AdmixCorr_',mark,' -e ',ex_param,' -i ',flockcommentDIR,'/SimDat',SDt,'/SimDatIn',SDt,' -o ',flockcommentDIR,'/SimDat',SDt,'/StructOuput_genos_slg_pipe.txt_dat00',SDt,'_k005_Rep00',R,'.txt > structureDumpola.txt', sep=""),wait=TRUE)     
    AdmixCorrRuntime[R]<-(proc.time()-ptm)[3]
    #NoAdmixCorr
    message("     NoAdmixCorr")
    ptm<-proc.time()
    system(paste(StructBinaryPath,' -m ',main_param,'NoAdmixCorr_',mark,' -e ',ex_param,' -i ',flockcommentDIR,'/SimDat',SDt,'/SimDatIn',SDt,' -o ',flockcommentDIR,'/SimDat',SDt,'/StructOuput_genos_slg_pipe.txt_dat00',SDt,'_k005_Rep01',R,'.txt  > structureDumpola.txt', sep=""),wait=TRUE)
    NoAdmixCorrRuntime[R]<-(proc.time()-ptm)[3]
    #NoAdmixNonCorr
    message("     NoAdmixNonCorr")
    ptm<-proc.time()
    system(paste(StructBinaryPath,' -m ',main_param,'NoAdmixNonCorr_',mark,' -e ',ex_param,' -i ',flockcommentDIR,'/SimDat',SDt,'/SimDatIn',SDt,' -o ',flockcommentDIR,'/SimDat',SDt,'/StructOuput_genos_slg_pipe.txt_dat00',SDt,'_k005_Rep02',R,'.txt  > structureDumpola.txt', sep=""),wait=TRUE)
    NoAdmixNonCorrRuntime[R]<-(proc.time()-ptm)[3]
  }
  write.table(x=AdmixCorrRuntime,file=paste(flockcommentDIR,'/SimDat',SDt,'/AdmixCorrRuntime.csv',sep=""),sep=",",quote=FALSE,append=TRUE,row.names=FALSE,col.names=FALSE)
  write.table(x=NoAdmixCorrRuntime,file=paste(flockcommentDIR,'/SimDat',SDt,'/NoAdmixCorrRuntime.csv',sep=""),sep=",",quote=FALSE,append=TRUE,row.names=FALSE,col.names=FALSE)
  write.table(x=NoAdmixNonCorrRuntime,file=paste(flockcommentDIR,'/SimDat',SDt,'/NoAdmixNonCorrRuntime.csv',sep=""),sep=",",quote=FALSE,append=TRUE,row.names=FALSE,col.names=FALSE)
}

##### PART 4: Run Distruct & Clumpp on datasets #### 
#we need to make an arena in ForStructure
system(paste('mkdir ',slg_pipeDIR,'/inputs/ForStructure/arena',sep=''))
#Arena - this is where we want to paste all our files lives
DistClumpFold<-paste(slg_pipeDIR,'/inputs/ForStructure/arena',sep="")

for (i in (1:(DatNum*Seedset))){
  #start by removing old directories 
  system(paste('rm ',DistClumpFold,'/StructOuput_genos_slg_pipe.txt* -print 2>/dev/null',sep=""))
  system(paste('rm -R ',slg_pipeDIR,'/inputs/ForStructure/clump_and_distruct/final_* -print 2>/dev/null',sep=""))
  system(paste('rm -R ',slg_pipeDIR,'/inputs/ForStructure/clump_and_distruct/intermediate -print 2>/dev/null',sep=""))
  
  #move each set of results for each Simulate dataset
  SysCmd<-paste('cp ',flockcommentDIR,'/SimDat',i,'/StructOuput_genos_slg_pipe.txt* ', DistClumpFold,sep="") #cp source destination
  system(SysCmd)
  
  #Then we just need to make sure a few files are correct
  write.table(Npops,file=paste(slg_pipeDIR,'/inputs/ForStructure/clump_and_distruct/num_pops.txt',sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
  write.table(paste(seq(from=1,to=Npops),seq(from=1,to=Npops),sep=" "),file=paste(slg_pipeDIR,'/inputs/ForStructure/clump_and_distruct/pop_names.txt',sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
  write.table(Nloci,file=paste(slg_pipeDIR,'/inputs/ForStructure/L.txt',sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
  write.table(N,file=paste(slg_pipeDIR,'/inputs/ForStructure/input/N.txt',sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
  
  #Now run clump and distruct!
  setwd(paste(slg_pipeDIR,'/inputs/ForStructure/clump_and_distruct',sep=""))
  ClmpDistCmd<- './script/ClumpAndDistructAll.sh 6'
  system(ClmpDistCmd)
  
  # move clump intermediate folder and all the ind pdfs
  cmd<-paste('cp -R ',slg_pipeDIR,'/inputs/ForStructure/clump_and_distruct/final_* ', flockcommentDIR,'/SimDat',i,sep="")
  system(cmd)
  cmd<-paste('cp -R ',slg_pipeDIR,'/inputs/ForStructure/clump_and_distruct/intermediate ', flockcommentDIR,'/SimDat',i,sep="")
  system(cmd)
}#over i simdat sets

##### PART 5: Compute Zero-One loss ####
index<-seq(from=1,to=((Reps*4)-(Reps-1)),by=Reps)

lossvc<-vector()
for (ds in (1:(DatNum*Seedset))){#how many data sets did we simulate
  for (r in (1:(Reps*4))){
    #first we need to assign each ind to a ref using its highest likelihood
    ff<-paste(flockcommentDIR,'/SimDat',ds,'/intermediate/Output_005.perms_',r,sep="")
    x<-read.table(ff)[c(6:10)]
    y<-apply(x,1,which.max)
    z<- cbind(read.table(ff)[4],y)
    STRUCTref<-unlist(z$y)
    y <- split(unlist(z$y), z$V4)
    
    ####Deal with label switiching - easiest way take the most common ref population given,
    #   but occasionally with really low Fst this won't work
    #   so on refs with the same most common label, we can use the one that we are more sure of
    #   ie. if label 3 is repeated and label 5 is missing which has a larger difference between 3 and 5
    
    ModeDist<-function(x){
      tl<- which.max(tabulate(x))
      return(tl)
    }
    
    Lab<-matrix(data=NA,nrow=5,ncol=2)
    Lab[,1]<-seq(1,5,1)
    Lab[,2]<-unlist(lapply(y,FUN=ModeDist))
    colnames(Lab)<-c('Known','STRUCTURE')
    
    # Are two references labelled the same?
    if(length(unique(Lab[,2]))==4){
      FC<-tabulate(Lab[,2],5)
      FCw<-which(FC!=1)# these are the ones that need to be assessed one might have >1 or 0
      FClrg<-which(FC>1)
      FClrgIndex<-unlist(lapply(FClrg,function(x) which(Lab[,2]==x)))#this is now a list
      FCsml<-which(FC<1)
      #which has the biggest difference between labels???
      FCb<-lapply(y[FClrgIndex],function(x) abs(diff((tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCw])))
      MSO<-FClrgIndex[which.max(FCb)]
      Lab[MSO,2]<-FClrg
      Lab[FClrgIndex[which.min(FCb)],2]<-FCsml
    }
    
    # if >2 sets of refs are labelled the same
    if((sum(tabulate(Lab[,2])==2))==3){
      FC<-tabulate(Lab[,2],5)
      FCw<-which(FC!=1)# these are the ones that need to be assessed one might have >1 or 0
      FClrg<-which(FC>1)
      FClrgIndex<-unlist(lapply(FClrg,function(x) which(Lab[,2]==x)))#this is now a list
      FCsml<-which(FC<1)
      #Edit for > 2 that are the same
      FCdifs<-c(FCsml[1],FClrg,FCsml[2])# this will be problematic if >3 diff. 
      FCb<-lapply(y[FClrgIndex],function(x) abs(diff((tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCdifs])))
      
      #which of the repeated reps has largest diff?
      MSO<-which.max(unlist(lapply(FCb,function(x) max(x))))
      #set that one
      Lab[MSO,2]<-FClrg # This just makes it stay the same...
      
      FClrg2index<-FClrgIndex[-MSO]
      # now assign the smaller ones and see if there is an issue
      MSO1<-lapply(y[FClrg2index],function(x) (tabulate(unlist(x))/sum(tabulate(unlist(x)))))
      MSO2<-lapply(y[FClrg2index],function(x) (tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCsml])
      MSOmax<-lapply(MSO2,function(x) max(x))
      MSOmax<-unlist(lapply(c(1,2),function(x) which(unlist(MSO1[x])==MSOmax[x])))
      Lab[FClrg2index,2]<-MSOmax
      
      if(length(unique(Lab[,2]))==4){
        #if both are assigned the same ref pop we then do what we did the first time!
        FC<-tabulate(Lab[,2],5)
        FCw<-which(FC!=1)# these are the ones that need to be assessed one might have >1 or 0
        FClrg<-which(FC>1)
        FClrgIndex<-which(Lab[,2]==FClrg)
        FCsml<-which(FC<1)
        
        #which has the biggest difference between labels???
        FCb<-lapply(y[FClrgIndex],function(x) abs(diff((tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCw])))
        MSO<-FClrgIndex[which.max(FCb)]
        Lab[MSO,2]<-FCw[(which.max(unlist(FCb)))]
        Lab[FClrgIndex[which.min(FCb)],3]<-FCw[(which.min(unlist(FCb)))]
      }#if
    }
    
    ## if 2 refs are assigned to 2 clusters
    if(length(unique(Lab[,2]))==3){
      FC<-tabulate(Lab[,2],5)
      FCw<-which(FC!=1)# these are the ones that need to be assessed one might have >1 or 0
      FClrg<-which(FC>1)
      FClrgIndex<-unlist(lapply(FClrg,function(x) which(Lab[,2]==x)))#this is now a list
      FCsml<-which(FC<1)
      
      FClrgIndex2<-lapply(FClrg,function(x) which(Lab[,2]==x))#this is now a list
      #assign both larges
      #1
      FCw<-c(FCsml[1],FClrg[1],FCsml[2]) # The first large and 2 unused
      FC1<-unlist(FClrgIndex2[1])
      FCb1<-lapply(y[FC1],function(x) mean(abs(diff((tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCw]))))
      FCb1Big<-which.max(FCb1)#this will give 1 or 2
      Lab[FC1[FCb1Big],2]<-FClrg[1]
      StillProb<-FC1[FC1!=FC1[FCb1Big]]
      #2
      FCw<-c(FCsml[1],FClrg[2],FCsml[2]) # The first large and 2 unused
      FC1<-unlist(FClrgIndex2[2])
      FCb1<-lapply(y[FC1],function(x) mean(abs(diff((tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCw]))))
      FCb1Big<-which.max(FCb1)#this will give 1 or 2
      Lab[FC1[FCb1Big],2]<-FClrg[2]
      StillProb2<-FC1[FC1!=FC1[FCb1Big]]
      
      SP<-c(StillProb,StillProb2)
      
      # well we should try to just assign to the larger for each...
      # SP needs to be assigned FCsml
      # just assign both to whichever is the largest and then if need be do the same thing as if we started over!
      Lab[SP[1],2]<-FCsml[which.max((tabulate(unlist(y[SP[1]]))/sum(tabulate(unlist(y[SP[1]]))))[FCsml])]
      Lab[SP[2],2]<-FCsml[which.max((tabulate(unlist(y[SP[2]]))/sum(tabulate(unlist(y[SP[2]]))))[FCsml])]
      
      # if there are three labelled the same thing... you can see that I found a lot of issues when
      # we decreased the population differentiation
      
      #if length lab ==4  
      if(length(unique(Lab[,2]))==4){
        #which has the biggest difference between labels???
        FC<-tabulate(Lab[,2],5)
        FCw<-which(FC!=1)# these are the ones that need to be assessed one might have >1 or 0
        FClrg<-which(FC>1)
        FClrgIndex<-which(Lab[,2]==FClrg)
        FCsml<-which(FC<1)
        FCb<-lapply(y[FClrgIndex],function(x) abs(diff((tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCw])))
        MSO<-FClrgIndex[which.max(unlist(FCb))]
        Lab[MSO,2]<-FClrg
        Lab[FClrgIndex[which.min(FCb)],2]<-FCsml
      }  
    }
    if(max(tabulate(Lab[,2],5))==3){
      #Edit for > 2 that are the same
      
      FCdifs<-c(FCsml[1],FClrg,FCsml[2])# this will be problematic if >3 diff. 
      FCb<-lapply(y[FClrgIndex],function(x) abs(diff((tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCdifs])))
      
      #which of the repeated reps has largest diff?
      MSO<-which.max(unlist(lapply(FCb,function(x) max(x))))
      #set that one
      Lab[MSO,2]<-FClrg # This just makes it stay the same...
      
      FClrg2index<-FClrgIndex[-MSO]
      # now assign the smaller ones and see if there is an issue
      MSO1<-lapply(y[FClrg2index],function(x) (tabulate(unlist(x))/sum(tabulate(unlist(x)))))
      MSO2<-lapply(y[FClrg2index],function(x) (tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCsml])
      MSOmax<-lapply(MSO2,function(x) max(x))
      MSOmax<-unlist(lapply(c(1,2),function(x) which(unlist(MSO1[x])==MSOmax[x])))
      Lab[FClrg2index,2]<-MSOmax
      if(length(unique(Lab[,2]))==2){
        #if both are assigned the same ref pop we then do what we did the first time!
        FC<-tabulate(Lab[,2],5)
        FCw<-which(FC!=1)# these are the ones that need to be assessed one might have >1 or 0
        FClrg<-which(FC>1)
        FClrgIndex<-which(Lab[,2]==FClrg)
        FCsml<-which(FC<1)
        
        #which has the biggest difference between labels???
        FCb<-lapply(y[FClrgIndex],function(x) abs(diff((tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCw])))
        MSO<-FClrgIndex[which.max(FCb)]
        Lab[MSO,2]<-FCw[(which.max(unlist(FCb)))]
        Lab[FClrgIndex[which.min(FCb)],2]<-FCw[(which.min(unlist(FCb)))]
      }#if
    }
    
    if(max(tabulate(Lab[,2]))>3){
      cat("I'm sorry, Dave. I'm afraid I can't do that. \n 
          You've created a dataset with an Fst that is so low that all programs have difficulties")
    }
    
    STRUCTrelab<-mapvalues(STRUCTref,Lab[,2],Lab[,1])
    
    #### Zero-One Loss
    # Let's compute the 0-1 loss for each run. We relabelled any mixups in the labeling above,
    # so it is straightforward using those labels and mach them to what exisits in each of the 
    # ref pops
    maxesFlockture<-Lab[,2]
    
    zero_one_loss <- sapply(1:5, function(h) {
      
      sum(maxesFlockture[h] != unlist(y[h]))
    })
    sum(zero_one_loss)
    lossvc[r]<-sum(zero_one_loss)
    }#over all r
  ZeroOneLosstab<-matrix(data=lossvc,byrow=FALSE,ncol=4,nrow=Reps)
  colnames(ZeroOneLosstab)<-c('A-C','NA-C','NA-NC','FLOCKTURE')
  rownames(ZeroOneLosstab)<-paste('Rep',seq(from=1,to=Reps,by=1))
  write.table(ZeroOneLosstab,paste(flockcommentDIR,"/SimDat",ds,"/ZeroOneLoss_SimDat",ds,".csv",sep=""),sep=",",col.names=TRUE,quote=FALSE)
# make and export boxplot for each model?

a<-melt(ZeroOneLosstab)
colnames(a)<-c('rep','model','loss')
pwFST<-as.character(round(read.table(paste(flockcommentDIR,'/SimDat',ds,'/AvgFst.txt',sep="")),4))
my.ylab = bquote(F[ST] == .(pwFST))
Ls<- ggplot(a, aes(x=model, y=loss)) + geom_boxplot() + labs(title=my.ylab)
  guides(fill=FALSE)+
  theme_bw() 


ggsave(plot=Ls,filename=paste(flockcommentDIR,'/SimDat',ds,'/LossCompSimDat',ds,'.pdf',sep=""),dpi=300,width=6,height=6,units='in')

#just see if there is a trend of qi values diverging as Fst increases
for (j in (1:Reps)){
structureNANCt<-read.table(paste(flockcommentDIR,'/SimDat',ds,'/intermediate','/Output_005.perms_',index[3]-1+j,sep=""))[,c(1,6:10)]
colnames(structureNANCt)<-c('ind',paste("r",seq(1,5,1),sep=""))
structureNANCt[,1]<-paste("Ind",structureNANCt[,1],sep="")
structureNANCt<-melt(structureNANCt)

flocktureT<-read.table(paste(flockcommentDIR,'/SimDat',ds,'/intermediate','/Output_005.perms_',index[4]-1+j,sep=""))[,c(1,6:10)]
colnames(flocktureT)<-c('ind',paste("r",seq(1,5,1),sep=""))
flocktureT[,1]<-paste("Ind",flocktureT[,1],sep="")
flocktureT<-melt(flocktureT)


#set up final dataset
Allq<-as.data.frame(cbind(flocktureT,structureNANCt[,3]))
colnames(Allq)<-c('ind','k','Flockture_q','Structure_NANC_q')

df<-as.data.frame(Allq)

p1<- ggplot(data=df, aes(y=Structure_NANC_q, x=Flockture_q, shape=factor(k))) +
  geom_point(alpha=1) +
  scale_shape_discrete(name="Cluster (k)") +
  scale_color_discrete(name="Cluster (k)") +
  #geom_rug(col=factor(k),alpha=.1)+
  theme_bw() +
  theme(legend.key = element_blank()) +
  theme(panel.grid=element_blank()) +
  theme(legend.position='bottom') +
  #theme(legend.background = element_rect(fill="gray95", size=.5, linetype="dotted")) +
  labs(x=expression(paste("FLOCKture ", italic(q[i])," values",sep=""), "values",sep=" "),y=expression(paste("STRUCTURE NA-NC ", italic(q[i])," values",sep=""))) +
  geom_abline(intercept=0, slope=1, linetype=2)
if(qi_indLoss==T){
ggsave(plot=p1,filename=paste(flockcommentDIR,'/SimDat',ds,'/QiCompRep',j,'SimDat',ds,'.pdf',sep=""),dpi=300,width=6,height=6,units='in')
}
}#over j reps
}#over all ds datasets

##### STEP 6: Make final figure ####

SumMat<-matrix(data=NA,ncol=3,nrow=Reps*4*DatNum)
for (i in 1:DatNum){
Fst<-round(unlist(rep(read.table(paste(flockcommentDIR,'/SimDat',i,'/AvgFst.txt',sep="")),Reps*4)),6)
LossMat<-read.csv(paste(flockcommentDIR,'/SimDat',i,'/ZeroOneLoss_SimDat',i,'.csv',sep=""),sep=",",colClasses=c('character',rep('numeric',4)))
LossMelt<-melt(LossMat)
SumMat[(((Reps*4*i)-(Reps*4)+1):(Reps*4*i)),(1:3)]<-as.matrix(cbind(LossMelt,Fst))
}
colnames(SumMat)<-c('Model','Loss','Fst')
SumMat<-as.data.frame(SumMat)
SumMat[,2]<-as.numeric(levels(SumMat[,2]))[SumMat[,2]]
SumMat<-SumMat[order(SumMat[,3]),]

PanLen<-nrow(SumMat)/4
SumMat2<-cbind(SumMat,rep(1:4,each=PanLen))
colnames(SumMat2)<-c('Model','Loss','Fst','Panel')
#need to rename levels
levels(SumMat2$Model)[5]<-'NA.C'
SumMat2[which(SumMat2[,1]=='A.NC'),1]<-'NA.C'

SumMat2$Model <- factor(SumMat2$Model, c("A.C", "NA.C", "NA.NC", "FLOCKTURE"))

P1<-ggplot(SumMat2,aes(x=Fst ,y=Loss,color=Model,fill=Model))  + 
  facet_wrap(~Panel,nrow=2,shrink=T,scales='free') + 
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white",outlier.colour = NA, 
               position = position_dodge(width=0.9))+
  theme(strip.background = element_blank(), strip.text = element_blank())+
  theme(axis.text.x = element_text(angle = 25, hjust = 1))+
  labs(x=expression('F'['ST']))

ggsave(plot=P1,filename=paste(flockcommentDIR,'/Loss',mark,'_Color.pdf',sep=''),dpi=300,width=11,height=8.5,units='in') 
ggsave(plot=P1+scale_colour_grey(start=.8,end=0),filename=paste(flockcommentDIR,'/Loss',mark,'_BW.pdf',sep=''),dpi=300,width=11,height=8.5,units='in')

#Plateau Table
setwd(flockcommentDIR)
Plata<-matrix(data=NA,ncol=Reps+1,nrow=(DatNum*Seedset))
#Lets Take a look at the Plateaus

for (d in 1:(DatNum*Seedset)){
  for (r in 1:Reps){
  Fst<-scan(file=paste(flockcommentDIR,'/SimDat',d,'/AvgFst.txt',sep=""))
  Plata[d,1]<-Fst
  plateau<-scan(file=paste(flockcommentDIR,'/SimDat',d,'/FlockturePlateaus.csv',sep=""),what="character",nlines=1,skip=(r-1))
  plateau<-paste(unlist(str_split(plateau,","))[which(unlist(str_split(plateau,","))>1)],collapse=",") 
  if(is.na(as.numeric(unlist(str_split(plateau,","))[1]))){plateau<-"1"}
  Plata[d,r+1]<-plateau
  }
}
colnames(Plata)<-c('Fst',paste('Rep',seq(1,Reps,1),sep=''))
Plata<-Plata[order(Plata[,1]),]
sink(file=paste('LatexPlateauTable_',mark,'.txt',sep=''),append=F)
print(xtable(Plata), include.rownames=FALSE)
sink()

system(paste("sed '/INSERTFILEHERE/r LatexPlateauTable_",mark,".txt' <Table2.tex >PlateuTable_",mark,"INT.tex",sep=""))
system(paste("sed 's/INSERTFILEHERE//' <PlateuTable_",mark,"INT.tex >PlateuTable_",mark,".tex",sep=""))
system("rm *INT.tex")            

##### clean up all the misc files sitting around
system("rm *DIR.txt") 
