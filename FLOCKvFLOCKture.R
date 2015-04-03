#### how does FLOCK compare with FLOCKTURE ####
# First we will need to set up some paths, so that we can 
# use the flockture script to analyze our datasets

flockcommentDIR <- normalizePath("../flock-compare") 
flocktureDIR <- normalizePath("../flockture") 

#Next, load some useful libraries
library(ggplot2)
library(stringr)
library(plyr)
library(reshape2)
library(miscFuncs)
library(xtable)
require("XLConnect")
options(java.parameters = "-Xmx4g" ) # if you have a lot of java heap space errors,
# you may need to run each of the Ds iterations individually, or run on a computer 
# than has more memory. 

# We will need to source the flockture code... it should be compiled first 
# by running 'make' in terminal within scr folder
setwd(flocktureDIR)
source("R/flockture.R")

#Now we will run Flockture#
# We ran FLOCK 9 times with 50 reps and 20 interations... make sure
# to system.time this to track how long it takes
Datasets <- 4 #how many datasets do we have to read in?
for (Ds in 1:Datasets){
  #Files are named 'BaseFile_SimX.txt' where X is the run number
  dat<-paste(flockcommentDIR,'/FLOCKvFLOCKTURE_data/data/BaseFile_Sim',Ds,'.txt',sep='')

# read in a data set:
D <- read.table(dat, row.names = 1)

# define some starting conditions if you want (this is useful to test initial conditions):
#sc <- c(c(rep(2,times=14),rep(1,times=15),rep(2,times=10),rep(1,times=11))
#)

FLCTR_PlateauRec<-vector()
Flockture.Runtime<-vector()
FlocktureReps<-9 #how many times are we running flockture
for (FR in 1:FlocktureReps){
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
Flockture.Runtime[FR]<-(proc.time()-ptm)[3]

# here is a vector of the different plateau lengths in sorted order:
sapply(psum$plateau_list, length)
FLCTR_PlateauRec1<-vector()
FLCTR_PlateauRec1<-sapply(psum$plateau_list, length)#[sapply(psum$plateau_list, length)>1]
FLCTR_PlateauRec[FR]<-paste(FLCTR_PlateauRec1,collapse=",")

#### Load some FLOCK RESULTS ####
FLOCKreps<-9
wd<-getwd()
Diff<-vector()
ZOlossFLOCK<-vector()
for (n in 1:FLOCKreps){

filename<-paste('Sim',Ds,'_k=5_run',n,'.xls',sep="")
excel.file <- file.path(flockcommentDIR,'FLOCKvFLOCKTURE_data','FlockResults',paste('Dat',Ds,sep=""),filename)
assign(x="FLOCKres",value=(readWorksheetFromFile(excel.file, sheet="All specs alloc and likelihoods", header=TRUE)))
FLOCKres<-FLOCKres[-(1),]#get rid of that pesky header that is imported
#FLOCK gives ref1 or ref2 in a column, lets pull ref out
Fk.ref<-as.numeric(gsub('ref','',FLOCKres[,2]))

# Which FLOCKTURE rep is the best?
# We should take the result with the longest plateau, 
# but if all plateaus are the same length we should take the one wiht the highest log.prob
# Longest plateau does not mean highest log.prob! So we need to do one and then the next.

#which reps belong to the largest plateau?
Max.pl<-max(psum$log_probs[psum$log_probs$iter==20,5])
Max.reps<-psum$log_probs[psum$log_probs$iter==20,][psum$log_probs[psum$log_probs$iter==20,5]==Max.pl,]
Best.rep<-Max.reps[which(Max.reps$log.prob==max(Max.reps$log.prob))[1],1]
  
Fkture.ref<-out$indivs[out$indivs$rep==Best.rep,4] 
#write.table(out$indivs[out$indivs$rep==Best.rep,],'FLOCKTUREqi.csv',sep=',',row.names=FALSE,col.names=T)

#### Comparison of FLOCK and FLOCKture ####
# Now we have to deal with some label switching.

# bind up Flock, Flockture and known value
Comp<-cbind(Fk.ref,Fkture.ref,c(rep(1,250),rep(2,250),rep(3,400),rep(4,550),rep(5,550)))
FLOCK_y <- split(Comp[,1], Comp[,3])
FLOCKture_y <-split(Comp[,2], Comp[,3])

# use apply function with a function to find the mode of each...
# we need this so we can use mapvalues to replace.... dealing with 
# label switching
ModeDist<-function(x){
 tl<- which.max(tabulate(x))
  return(tl)
}

Lab<-matrix(data=NA,nrow=5,ncol=3)
Lab[,1]<-seq(1,5,1)
Lab[,2]<-unlist(lapply(FLOCK_y,FUN=ModeDist))
Lab[,3]<-unlist(lapply(FLOCKture_y,FUN=ModeDist))
colnames(Lab)<-c('Known','FLOCK','FLOCKTURE')

# what if the program does really poorly and the mode isn't useful
# for determining which ref group is which?

#for flock
FC<-tabulate(Lab[,2],5)
FCw<-which(FC!=1)# these are the ones that need to be assessed one might have >1 or 0
FClrg<-which(FC>1)
FClrgIndex<-which(Lab[,2]==FClrg)
FCsml<-which(FC<1)

#if Only 2 refs are assigned the same number
if(length(unique(Lab[,2]))==4){ 

#which has the biggest difference between labels???
FCb<-lapply(FLOCK_y[FClrgIndex],function(x) abs(diff((tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCw])))
MSO<-FClrgIndex[which.max(FCb)]
Lab[MSO,2]<-FCw[(which.max(unlist(FCb)))]
Lab[FClrgIndex[which.min(FCb)],2]<-FCw[(which.min(unlist(FCb)))]
}


# if 3 are three that are assigned the same?
if(max(tabulate(Lab[,2],5))==3){
  FCdifs<-c(FCsml[1],FClrg,FCsml[2])# this will be problematic if >3 diff. 
  FCb<-lapply(FLOCKture_y[FClrgIndex],function(x) abs(diff((tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCdifs])))
  
  #which of the repeated reps has largest diff?
  MSO<-which.max(unlist(lapply(FCb,function(x) max(x))))
  #set that one
  Lab[MSO,2]<-FClrg # This just makes it stay the same...
  
  FClrg2index<-FClrgIndex[-MSO]
  # now assign the smaller ones and see if there is an issue
  MSO1<-lapply(FLOCKture_y[FClrg2index],function(x) (tabulate(unlist(x))/sum(tabulate(unlist(x)))))
  MSO2<-lapply(FLOCKture_y[FClrg2index],function(x) (tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCsml])
  MSOmax<-lapply(MSO2,function(x) max(x))
  MSOmax<-lapply(c(1,2),function(x) which(unlist(MSO1[x])==MSOmax[x]))
  Lab[FClrg2index,2]<-c(MSOmax[[1]][1],MSOmax[[2]][1])
  if(max(tabulate(Lab[,2]))==2){
    #if both are assigned the same ref pop we then do what we did the first time!
    FC<-tabulate(Lab[,2],5)
    FCw<-which(FC!=1)# these are the ones that need to be assessed one might have >1 or 0
    FClrg<-which(FC>1)
    FClrgIndex<-which(Lab[,2]==FClrg)
    FCsml<-which(FC<1)
    
    #which has the biggest difference between labels???
    FCb<-lapply(FLOCKture_y[FClrgIndex],function(x) abs(diff((tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCw])))
    MSO<-FClrgIndex[which.max(FCb)]
    Lab[MSO,2]<-FCw[(which.max(unlist(FCb)))]
    Lab[FClrgIndex[which.min(FCb)],2]<-FCw[(which.min(unlist(FCb)))]
  }#if
}#if
######
if(max(tabulate(Lab[,2]))>3){
  cat("Analysis yielded fewer than five clusters.")
}



################################################################################################
#for flockture
#if(length(unique(Lab[,3]))<5){
  #sort(as.numeric(unique(Lab[,2])))# does flock have all 5 group labels?
FC<-tabulate(Lab[,3],5)
FCw<-which(FC!=1)# these are the ones that need to be assessed one might have >1 or 0
FClrg<-which(FC>1)
FClrgIndex<-unlist(lapply(FClrg,function(x) which(Lab[,3]==x)))#this is now a list
FCsml<-which(FC<1)
# Are two references labelled the same?
if(length(unique(Lab[,3]))==4){
  #which has the biggest difference between labels???
  FCb<-lapply(FLOCKture_y[FClrgIndex],function(x) abs(diff((tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCw])))
  MSO<-FClrgIndex[which.max(FCb)]
  Lab[MSO,3]<-FClrg
  Lab[FClrgIndex[which.min(FCb)],3]<-FCsml
}
  
# if >2 sets of refs are labelled the same
if((sum(tabulate(Lab[,3])==2))==3){
  #Edit for > 2 that are the same
  FCdifs<-c(FCsml[1],FClrg,FCsml[2])# this will be problematic if >3 diff. 
  FCb<-lapply(FLOCKture_y[FClrgIndex],function(x) abs(diff((tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCdifs])))
  
  #which of the repeated reps has largest diff?
  MSO<-which.max(unlist(lapply(FCb,function(x) max(x))))
  #set that one
  Lab[MSO,3]<-FClrg # This just makes it stay the same...
  
  FClrg2index<-FClrgIndex[-MSO]
  # now assign the smaller ones and see if there is an issue
  MSO1<-lapply(FLOCKture_y[FClrg2index],function(x) (tabulate(unlist(x))/sum(tabulate(unlist(x)))))
  MSO2<-lapply(FLOCKture_y[FClrg2index],function(x) (tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCsml])
  MSOmax<-lapply(MSO2,function(x) max(x))
  MSOmax<-unlist(lapply(c(1,2),function(x) which(unlist(MSO1[x])==MSOmax[x])))
  Lab[FClrg2index,3]<-MSOmax
   
  if(length(unique(Lab[,3]))==4){
#if both are assigned the same ref pop we then do what we did the first time!
FC<-tabulate(Lab[,3],5)
FCw<-which(FC!=1)# these are the ones that need to be assessed one might have >1 or 0
FClrg<-which(FC>1)
FClrgIndex<-which(Lab[,3]==FClrg)
FCsml<-which(FC<1)

#which has the biggest difference between labels???
FCb<-lapply(FLOCKture_y[FClrgIndex],function(x) abs(diff((tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCw])))
MSO<-FClrgIndex[which.max(FCb)]
Lab[MSO,3]<-FCw[(which.max(unlist(FCb)))]
Lab[FClrgIndex[which.min(FCb)],3]<-FCw[(which.min(unlist(FCb)))]
  }#if
}

###### if 2 refs are assigned to 2 clusters
if(length(unique(Lab[,3]))==3){
  FClrgIndex2<-lapply(FClrg,function(x) which(Lab[,3]==x))#this is now a list
  #assign both larges
  #1
  FCw<-c(FCsml[1],FClrg[1],FCsml[2]) # The first large and 2 unused
  FC1<-unlist(FClrgIndex2[1])
  FCb1<-lapply(FLOCKture_y[FC1],function(x) mean(abs(diff((tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCw]))))
  FCb1Big<-which.max(FCb1)#this will give 1 or 2
  Lab[FC1[FCb1Big],3]<-FClrg[1]
  StillProb<-FC1[FC1!=FC1[FCb1Big]]
  #2
  FCw<-c(FCsml[1],FClrg[2],FCsml[2]) # The first large and 2 unused
  FC1<-unlist(FClrgIndex2[2])
  FCb1<-lapply(FLOCKture_y[FC1],function(x) mean(abs(diff((tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCw]))))
  FCb1Big<-which.max(FCb1)#this will give 1 or 2
  Lab[FC1[FCb1Big],3]<-FClrg[2]
  StillProb2<-FC1[FC1!=FC1[FCb1Big]]
  
  SP<-c(StillProb,StillProb2)
  
  # well we should try to just assign to the larger for each...
  # SP needs to be assigned FCsml
  # just assign both to whichever is the largest and then if need be do the same thing as if we started over!
  Lab[SP[1],3]<-FCsml[which.max((tabulate(unlist(FLOCKture_y[SP[1]]))/sum(tabulate(unlist(FLOCKture_y[SP[1]]))))[FCsml])]
  Lab[SP[2],3]<-FCsml[which.max((tabulate(unlist(FLOCKture_y[SP[2]]))/sum(tabulate(unlist(FLOCKture_y[SP[2]]))))[FCsml])]

#if there are three labelled the same thing

#if length lab ==4  
if(length(unique(Lab[,3]))==4){
  #which has the biggest difference between labels???
  FC<-tabulate(Lab[,3],5)
  FCw<-which(FC!=1)# these are the ones that need to be assessed one might have >1 or 0
  FClrg<-which(FC>1)
  FClrgIndex<-which(Lab[,3]==FClrg)
  FCsml<-which(FC<1)
  FCb<-lapply(FLOCKture_y[FClrgIndex],function(x) abs(diff((tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCw])))
  MSO<-FClrgIndex[which.max(unlist(FCb))]
  Lab[MSO,3]<-FClrg
  Lab[FClrgIndex[which.min(FCb)],3]<-FCsml
}  
}
if(max(tabulate(Lab[,3],5))==3){
  #Edit for > 2 that are the same
  
  FCdifs<-c(FCsml[1],FClrg,FCsml[2])# this will be problematic if >3 diff. 
  FCb<-lapply(FLOCKture_y[FClrgIndex],function(x) abs(diff((tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCdifs])))
  
  #which of the repeated reps has largest diff?
  MSO<-which.max(unlist(lapply(FCb,function(x) max(x))))
  #set that one
  Lab[MSO,3]<-FClrg # This just makes it stay the same...
  
  FClrg2index<-FClrgIndex[-MSO]
  # now assign the smaller ones and see if there is an issue
  MSO1<-lapply(FLOCKture_y[FClrg2index],function(x) (tabulate(unlist(x))/sum(tabulate(unlist(x)))))
  MSO2<-lapply(FLOCKture_y[FClrg2index],function(x) (tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCsml])
  MSOmax<-lapply(MSO2,function(x) max(x))
  MSOmax<-unlist(lapply(c(1,2),function(x) which(unlist(MSO1[x])==MSOmax[x])))
  Lab[FClrg2index,3]<-MSOmax
  if(length(unique(Lab[,3]))==2){
    #if both are assigned the same ref pop we then do what we did the first time!
    FC<-tabulate(Lab[,3],5)
    FCw<-which(FC!=1)# these are the ones that need to be assessed one might have >1 or 0
    FClrg<-which(FC>1)
    FClrgIndex<-which(Lab[,3]==FClrg)
    FCsml<-which(FC<1)
    
    #which has the biggest difference between labels???
    FCb<-lapply(FLOCKture_y[FClrgIndex],function(x) abs(diff((tabulate(unlist(x))/sum(tabulate(unlist(x))))[FCw])))
    MSO<-FClrgIndex[which.max(FCb)]
    Lab[MSO,3]<-FCw[(which.max(unlist(FCb)))]
    Lab[FClrgIndex[which.min(FCb)],3]<-FCw[(which.min(unlist(FCb)))]
  }#if
}

if(length(unique(Lab[,3]))==2){
  cat("Only two reference groups!!!")
}
#############################################################################################
FLOCKrelab<-mapvalues(Fk.ref,Lab[,2],Lab[,1])
FLOCKturerelab<-mapvalues(Fkture.ref,Lab[,3],Lab[,1])

#how many differences?
Diff[n]<-length(FLOCKrelab)-sum(FLOCKrelab==FLOCKturerelab)

####Zero-One Loss####

# we only need to do this once for each dataset
maxesFLOCK <- Lab[,2]
maxesFlockture<-Lab[,3]

zero_one_loss <- sapply(1:5, function(h) {
  c(sum(maxesFLOCK[h] != unlist(FLOCK_y[h])),
  sum(maxesFlockture[h] != unlist(FLOCKture_y[h])))
})
if(FR==1){
write.table(x=sum(zero_one_loss[1,]),file='ZOlossFLOCK.txt',append=TRUE,row.names=FALSE,col.names=FALSE)
} #if FR==1
if(n==1){ # we only need to write out flockture Z-O loss when we change F
  write.table(x=sum(zero_one_loss[2,]),file='ZOlossFLOCKture.txt',sep=",",quote=FALSE,append=TRUE,row.names=FALSE,col.names=FALSE)
} # if n==1
####Plot all of the Qi values for each of the 5 clusters####
# if they behave similarly then they will be lined up on the 1:1 line
# data should be in ind, k, flock_q, flockture_q columns

FLOCKTUREtidy<-as.data.frame(matrix(data=NA,ncol=3,nrow=2000*5,))
FLOCKTUREtidy[,1]<-as.numeric(unlist(out$indivs[out$indivs$rep==Best.rep,5:9]))
FLOCKTUREtidy[,2]<-rep(1:5,each=2000)
FLOCKTUREtidy[,3]<-rep(FLOCKres[,1],5)
FLOCKTUREtidy[,2]<-mapvalues(FLOCKTUREtidy[,2],Lab[,3],Lab[,1])#relabel the columns
FLOCKTUREts<-FLOCKTUREtidy[ order(FLOCKTUREtidy[,2]), ]

FLOCKtidy<-as.data.frame(matrix(data=NA,ncol=3,nrow=2000*5))
FLOCKtidy[,1]<-as.numeric(unlist(FLOCKres[,(15:19)]))
FLOCKtidy[,2]<-rep(1:5,each=2000)
FLOCKtidy[,3]<-rep(FLOCKres[,1],5)
FLOCKtidy[,2]<-mapvalues(FLOCKtidy[,2],Lab[,2],Lab[,1])#relabel the columns
FLOCKts<-FLOCKtidy[ order(FLOCKtidy[,2]), ]


#set up final dataset
Allq<-as.data.frame(cbind(FLOCKTUREts[,3],FLOCKTUREts[,1],FLOCKts[,1:2]))
colnames(Allq)<-c('ind','Flockture_q','Flock_q','k')

df<-as.data.frame(Allq)

p1 <- ggplot(data=df, aes(x=Flock_q, y=Flockture_q, shape=factor(k))) +
  geom_point(alpha=0.5) +
  scale_shape_discrete(name="Cluster (k)") +
  scale_color_discrete(name="Cluster (k)") +
  #geom_rug(col=factor(k),alpha=.1)+
  theme_bw() +
  theme(legend.key = element_blank()) +
  theme(panel.grid=element_blank()) +
  theme(legend.position='bottom') +
  #theme(legend.background = element_rect(fill="gray95", size=.5, linetype="dotted")) +
  labs(y=expression(paste("FLOCKture ", italic(q[i])," values",sep=""), "values",sep=" "),x=expression(paste("FLOCK ", italic(q[i])," values",sep=""))) +
  geom_abline(intercept=0, slope=1, linetype=2)

ggsave(plot=p1,filename=paste('FigOut_FLCKTR',FR,'_FLOCK',n,'.pdf',sep=""),dpi=300,width=6,height=6,units='in')

#### FLOCK plateau lengths #####
# These are stored in a sep file from all the other information
# so we will have to pull them in one at a time.
#only need to do this once too!
if(FR==1){
filename2<-paste('Sim',Ds,'_run',n,'_plateaurecord.xls',sep="")
excel.file2 <- file.path(flockcommentDIR,'FLOCKvFLOCKTURE_data','FlockResults',paste('Dat',Ds,sep=""),filename2)
assign(x="FLOCKpl",value=(readWorksheetFromFile(excel.file2, sheet="Sheet1", header=TRUE)))
write.table(x=FLOCKpl[,4],file='FlockPlateaus.csv',sep=",",quote=FALSE,append=TRUE,row.names=FALSE,col.names=FALSE)
}# if - write out the plateau lengths to file as they are read in. 
}#over n reps of flock
Diff.avg<-paste(round(mean(Diff),2),' (',round(sd(Diff),2),')',sep="",collapse="")
write.table(x=Diff.avg,file='AvgDifferences.csv',sep=",",quote=FALSE,append=TRUE,row.names=FALSE,col.names=FALSE)
}#over FR reps of flockture
write.table(x=FLCTR_PlateauRec,file='FlockturePlateaus.csv',sep=",",quote=FALSE,append=TRUE,row.names=FALSE,col.names=FALSE)
write.table(x=Flockture.Runtime,file='FlocktureRuntime.csv',sep=",",quote=FALSE,append=TRUE,row.names=FALSE,col.names=FALSE)
####

####Let's make a table #####
tab<-matrix(data=NA,nrow=9,ncol=8)
tab[,1]<-seq(1,9,1)
FRTfile<-file.path(flockcommentDIR,'FLOCKvFLOCKTURE_data','FlockResults',paste('Dat',Ds,sep=""),paste('FlockRuntimesDat',Ds,'.csv',sep=''))
tab[,2]<-as.numeric(gsub('0:','',read.table(FRTfile,sep=",")[,3]))*60
tab[,3]<-gsub(" ", "", readLines('FlockPlateaus.csv'))
tab[,4]<-read.table('ZOlossFLOCK.txt',sep=",")[,1]
tab[,5]<-round((read.table('FlocktureRuntime.csv',sep=",")[,1]),2)
if(as.numeric(str_split(readLines('FlockturePlateaus.csv'),pattern=",")[[1]][1])>1){
  tab[,6]<-gsub(',{2,}','',gsub('1','',readLines('FlockturePlateaus.csv')))               
}else{tab[,6]<-rep(0,9)}
tab[,7]<-read.table('ZOlossFLOCKture.txt',sep=",")[,1]
tab[,8]<-as.character(read.table('AvgDifferences.csv',sep=",")[,1])
write.table(tab,file=paste('FLOCKvFLOCKTUREdat',Ds,'.csv',sep=''),sep=',',row.names=FALSE,col.names=FALSE)

#### Move all the files we make into the right place
Dest<-file.path(flockcommentDIR,'FLOCKvFLOCKTURE_data','FlockResults',paste('Dat',Ds,sep=''),'Comp')
#make a directory
system(command=paste('mkdir ',Dest,sep=''))
system(command= paste('mv ',paste(wd,'/FigOut*',sep=''),paste(Dest)))
system(command= paste('mv ',paste(wd,'/*Plateaus.csv',sep=''),paste(Dest)))
system(command= paste('mv ',paste(wd,'/ZOloss*',sep=''),paste(Dest)))
system(command= paste('mv ',paste(wd,'/FlocktureRuntime*',sep=''),paste(Dest)))
system(command= paste('mv ',paste(wd,'/AvgDiff*',sep=''),paste(Dest)))
system(command= paste('mv ',paste(wd,'/FLOCKvFLOCKTURE*',sep=''),paste(Dest)))
cat(paste('Dataset ',Ds,' done and results located in Comp folder. \n',sep=''))
}#over Ds



#### Table to compare FLOCK and FLOCKTURE ####
#read in all the results
ext<-sapply(1:4, function(x) paste('read.csv("',flockcommentDIR,'/FLOCKvFLOCKTURE_data/FlockResults/Dat',x,'/Comp/FLOCKvFLOCKTUREdat',x,'.csv",','na.strings="NA",header=FALSE)',sep=""))
table1<-rbind(eval(parse(text=ext[1])),eval(parse(text=ext[2])),eval(parse(text=ext[3])),eval(parse(text=ext[4])))
table1[]<-lapply(table1,as.character)
table1[is.na(table1)]<-0
table1<-cbind(rep(1:4,each=9),table1)
colnames(table1)<-c('Sim','Rep','Run Time','Plateau Record','Zero-One Loss','Run Time','Plateau Record','Zero-One Loss','Avg. Differences (SD)')

setwd(paste(flockcommentDIR,'/FLOCKvFLOCKTURE_data',sep=""))
sink('Table1_content.tex')
print(xtable(table1), include.rownames=FALSE)
sink()
sink()
#Use sed to paste into a better fomatted Latex document!
#sed '/INCLUDE/ r file' <in >out 
##sed '/INCLUDE/d' file >out 

#which file to insert
    sapply(1:36,function(x) write.table(unlist(scan('Table1_content.tex',skip=7+x,nlines=1,what='character',sep='\n')),file='Nums.tex',append=T,row.names=FALSE,col.names=FALSE,quote=FALSE))
    path<-"sed '/INSERT NUMBERS HERE/ r Nums.tex' <~/Flockture/FlockResults/Table1Skeleton.tex >Intermediate.tex"
    system(path)
    system("sed '/INSERT NUMBERS HERE/d' 'Intermediate.tex' >Tabel1.tex")
    
#clean up all the .tex files not in use
system('rm Intermediate.tex')
system('rm Table1_content.tex')
system('rm Nums.tex')

#### Create Loss Figure ####
facet<-rep(rep(c('a','b','c','d'),each=9),2)
c<-cbind(c,facet)
theme_set(theme_gray(base_size = 18))
LossFvFl<-ggplot(c,aes(x=Model ,y=as.numeric(ModelLoss)))  + 
  facet_wrap(~facet,ncol=2,scales='free') + 
  geom_boxplot(fill="white",outlier.colour = NA)+
  theme(axis.title.x = element_text(vjust=-0.75,size=20),axis.text=element_text(size=14,colour="black")) +
  labs(x='Model',y='Loss')
ggsave(plot=LossFvFl,filename='LossFlockvFlockture.pdf',dpi=300,width=8,height=8,units='in') 
