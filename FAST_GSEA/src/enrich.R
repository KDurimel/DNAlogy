######################################################################################################################
############# GO-TERMS ENRICHMENT
######################################################################################################################

#######################################################
#GENERAL SETTINGS
#######################################################

#Load GOdb silently
suppressMessages(library(GO.db))

# Raise number of printable strings in order to be able to capture large outputs
options(max.print=99999999)

# parse arguments from master script
args <- commandArgs(TRUE)

outputPrefix=args[1]    
outputPrefix=gsub(" ", "",outputPrefix, fixed = TRUE)  
ech_tmp=args[2]
univ_tmp=args[3]
# read GO mapping files from master script
go_ech_file=readLines(ech_tmp)
go_univ_file=readLines(univ_tmp)



#####################################################################
# BROWSE THE GO-TEMRS ONTOLOGY TREE
#####################################################################

################### ASPECT 1: CC ################################

###RETRIEVE UNIVERSE GO-TERMS ANCESTORS###

print('Browsing sample ontology tree for cellular components...')

capture.output(for(i in 1:length(go_univ_file))
{
  aspect=1
  #Try, in order to avoid errors if GO-term is not a CC or dont have any ancestors
  try(print(aspect<-get(go_univ_file[i],GOCCANCESTOR)),TRUE)
  options(warn=-1) #This test is valid but generate msgerr, so, turn-off msgerr
  if(aspect!=1) #if aspect!=1 aspect is an GOCC
  {
    print(go_univ_file[i]) #print first GO-term if is its a CC
  }
  options(warn=0) # test ok, turn-on msgerr
},file=(paste(outputPrefix,"/gocc_univ.txt",sep="")))

#Open file and keep GO terms only (update file and erase old version)
data=scan(paste(outputPrefix,"/","gocc_univ.txt",sep=""),what="character"())#opened with fonction scan (+as character) in order to grep efficiently with go_only

#Special function in order to dont show columns and row names during the writing of the output file
print.go_only <- function(m)
{
  write.table(format(m, justify="right"),
              row.names=F, col.names=F, quote=F)
}

#Save results which only contains GO Id's
capture.output(print.go_only(grep(pattern ="GO:" , data, value = TRUE, fixed = TRUE)),file=paste(outputPrefix,"/gocc_univ.txt",sep=""))

###RETRIEVE SAMPLE GO-TERMS ANCESTORS###
print('Browsing universe ontology tree for cellular components...')
capture.output(for(i in 1:length(go_ech_file))
{
  #Try, in order to avoid errors if GO-term is not a CC or dont have any ancestors
  aspect=1
  #Try, in order to avoid errors if GO-term is not a CC or dont have any ancestors
  try(print(aspect<-get(go_ech_file[i],GOCCANCESTOR)),TRUE)
  options(warn=-1)
  if(aspect!=1) #if aspect!=1 aspect is an GOCC
  {
    print(go_ech_file[i]) #print first GO-term if is its a CC
  }
  options(warn=0)
},file=(paste(outputPrefix,"/gocc_ech.txt",sep="")))

#Open file and keep GO terms only (update file and erase old version)
data=scan(paste(outputPrefix,"/gocc_ech.txt",sep=""),what="character"())#opened with fonction scan (+as character) in order to grep efficiently with go_only

#Save results which only contains GO Id's
capture.output(print.go_only(grep(pattern ="GO:" , data, value = TRUE, fixed = TRUE)),file=(paste(outputPrefix,"/gocc_ech.txt",sep="")))

################### ASPECT 2: BP ################################

###RETRIEVE UNIVERSE GO-TERMS ANCESTORS###

print('Browsing universe ontology tree for biological processes...')

capture.output(for(i in 1:length(go_univ_file))
{
  aspect=1
  #Try, in order to avoid errors if GO-term is not a CC or dont have any ancestors
  try(print(aspect<-get(go_univ_file[i],GOBPANCESTOR)),TRUE)
  options(warn=-1)
  if(aspect!=1) #if aspect!=1 aspect is an GOCC
  {
    print(go_univ_file[i]) #print first GO-term if is its a CC
  }
  options(warn=0)
},file=(paste(outputPrefix,"/gobp_univ.txt",sep="")))

#Open file and keep GO terms only (update file and erase old version)
data=scan(paste(outputPrefix,"/gobp_univ.txt",sep=""),what="character"()) 
capture.output(print.go_only(grep(pattern = "GO:" , data, value = TRUE, fixed = TRUE)),file=(paste(outputPrefix,"/gobp_univ.txt",sep="")))

###RETRIEVE SAMPLE GO-TERMS ANCESTORS###
print('Browsing sample ontology tree for biological processes...')
capture.output(for(i in 1:length(go_ech_file))
{
  #Try, in order to avoid errors if GO-term is not a CC or dont have any ancestors
  aspect=1
  #Try, in order to avoid errors if GO-term is not a CC or dont have any ancestors
  try(print(aspect<-get(go_ech_file[i],GOBPANCESTOR)),TRUE)
  options(warn=-1)
  if(aspect!=1) #if aspect!=1 aspect is an GOCC
  {
    print(go_ech_file[i]) #print first GO-term if is its a CC
  }
  options(warn=0)
},file=(paste(outputPrefix,"/gobp_ech.txt",sep="")))

#Open file and keep GO terms only (update file and erase old version)
data=scan(paste(outputPrefix,"/gobp_ech.txt",sep=""),what="character"())#opened with fonction scan (+as character) in order to grep efficiently with go_only
#Save results which only contains GO Id's
capture.output(print.go_only(grep(pattern ="GO:" , data, value = TRUE, fixed = TRUE)),file=(paste(outputPrefix,"/gobp_ech.txt",sep="")))


################### ASPECT 3: MF ################################

###RETRIEVE UNIVERSE GO-TERMS ANCESTORS###

print('Browsing universe ontology tree for molecular functions...')

capture.output(for(i in 1:length(go_univ_file))
{
  #Try, in order to avoid errors if GOterm is not a MF or dont have any ancestors
  aspect=1
  #Try, in order to avoid errors if GO-term is not a CC or dont have any ancestors
  try(print(aspect<-get(go_univ_file[i],GOMFANCESTOR)),TRUE)
  options(warn=-1)
  if(aspect!=1) #if aspect!=1 aspect is an GOCC
  {
    print(go_univ_file[i]) #print first GO-term if is its a CC
  }
  options(warn=0)
},file=(paste(outputPrefix,"/gomf_univ.txt",sep="")))

#Open file and keep GO terms only (update file and erase old version)
data=scan(paste(outputPrefix,"/gomf_univ.txt",sep=""),what="character"())
capture.output(print.go_only(grep(pattern = "GO:" , data, value = TRUE, fixed = TRUE)),file=(paste(outputPrefix,"/gomf_univ.txt",sep="")))

###RETRIEVE SAMPLE GO-TERMS ANCESTORS###
print('Browsing sample ontology tree for molecular functions...')
capture.output(for(i in 1:length(go_ech_file))
{
  #Try, in order to avoid errors if GO-term is not a CC or dont have any ancestors
  aspect=1
  #Try, in order to avoid errors if GO-term is not a CC or dont have any ancestors
  try(print(aspect<-get(go_ech_file[i],GOMFANCESTOR)),TRUE)
  options(warn=-1)
  if(aspect!=1) #if aspect!=1 aspect is an GOCC
  {
    print(go_ech_file[i]) #print first GO-term if is its a CC
  }
  options(warn=0)
},file=(paste(outputPrefix,"/gomf_ech.txt",sep="")))

#Open file and keep GO terms only (update file and erase old version)
data=scan(paste(outputPrefix,"/gomf_ech.txt",sep=""),what="character"())#opened with fonction scan (+as character) in order to grep efficiently with go_only
#Save results which only contains GO Id's
capture.output(print.go_only(grep(pattern ="GO:" , data, value = TRUE, fixed = TRUE)),file=(paste(outputPrefix,"/gomf_ech.txt",sep="")))


######################################################################################################################
############# GO-THREE WALKING AND RELATED STATISTICAL TESTS 
######################################################################################################################

print('Performing hypergeometric tests...')

############################################
#HYPERGEOMETRIC TEST
############################################

#phyper(a,b,c,d,lower.tail = TRUE)
#
#a : GO-term hits in sample
#b : GO-term hits in universe
#c : universe of GO-terms length - GO-terms sample length
#d : GO-terms sample length
#H0: overrepresented GO-term
#H1: not overrepresented GO-term
#Lower.tail = TRUE : invert HO and H1 hypothesis.

#########################################################
#INITIALIZE VALUES FOR UNIVERSE AND DATASETS 
#########################################################
univers_mf=read.table(paste(outputPrefix,"/gomf_univ.txt",sep=""))
echantillon_mf=read.table(paste(outputPrefix,"/gomf_ech.txt",sep=""))

univers_cc=read.table(paste(outputPrefix,"/gocc_univ.txt",sep=""))
echantillon_cc=read.table(paste(outputPrefix,"/gocc_ech.txt",sep=""))

univers_bp=read.table(paste(outputPrefix,"/gobp_univ.txt",sep=""))
echantillon_bp=read.table(paste(outputPrefix,"/gobp_ech.txt",sep=""))
#######################################################

#Special function in order to do not keep undesired characters
onlythis<- function(m)
{
  write.table(format(m, justify="right"),
              row.names=F, col.names=F, quote=F)
}

########################################################
#R SETTINGS AND SOME VALUES SETTINGS
########################################################

Go_univ=univers_mf #Universe
Go_ech=echantillon_mf #Node
options(max.print=99999999) #Raise print capacity
Go_uniq=Go_univ[!duplicated(Go_univ),] #Go_uniq used for loop in order to dont retrieve a goterm multiple times if its duplicated
parent="all" #Root node can be considered as "all" in Go.db

#########################################################
#COMPUTING AND WRITE RESULTS IN FILE FOR EACH ASPECT:
#########################################################

#################### MOLECULAR FUNCTION #################

#####PHYPER TEST WITH BONFERRONI CORRECTION#####
pvals=""
for(i in 1:length(Go_uniq))
{
  targeted_goterm=toString(Go_uniq[i]) #tostring in order to request value
  pvals[i]<-phyper(length(which(Go_ech==targeted_goterm)),length(which(Go_univ==targeted_goterm)),length(Go_univ[,1])-length(which(Go_univ==targeted_goterm)),length(Go_ech[,1]),lower.tail = FALSE) 
}
pvals2<-as.vector(pvals)
pvals_adjusted<-p.adjust(pvals2, method="bonferroni", n=length(pvals))
a=1

######PHYPER OK, NOW ENHANCE PHYPER RESULTS#####
capture.output(
  for(i in 1:length(Go_uniq))
  {
    targeted_goterm=toString(Go_uniq[i]) #tostring in order to request value
    #Targeted_goterm: GOterm currently processed##Number of hits##Expected number of hits##phyper(count of this GO in the batch,count of this go in the universe, universe length,batch length)
    expectednumberofhits=(((length(which(Go_univ==targeted_goterm)))/(length(Go_univ[,1])))*length(Go_ech[,1]))
    trm=Term(targeted_goterm)
    ##########################################################################
    #RETRIEVE GO-TERM LEVEL : by parents "walking" method (is-a , part-of)
    ##########################################################################
    #First request
    parents=get(targeted_goterm,GOMFPARENTS) 
    #If parents founded-->next DAG level exists.So repeat request with the parents of the next DAG level
    if(parents[1]!=parent)
    {
      level=1 #level in the DAG
      tryCatch({
        for(i in 1:1000) #Arbitrary.
        {
          parents=get(parents[1],GOMFPARENTS) #[1] because prioritize for "is-a"
          level=level+1
          if(parents[1]==parent) stop()#If node ="all" , all go three processed,no more DAG levels.Stop the requests loop
        }
      }, error=function(e){cat()})#Nothing in message error in order to dont print msgerros in results file
    }
    #If no parents founded-->term is on highest level. This is Aspect term-->lvl
    else
    {
      level=1
    }
    #Print results foreach go-term
    if(length(which(Go_ech==targeted_goterm))!=0) #if go-term is present in sample
    {
	    col3=as.character(length(which(Go_ech==targeted_goterm)))
	    col6=as.character(phyper(length(which(Go_ech==targeted_goterm)),length(which(Go_univ==targeted_goterm)),length(Go_univ[,1])-length(which(Go_univ==targeted_goterm)),length(Go_ech[,1]),lower.tail = FALSE))
	    cat(targeted_goterm,trm,col3,expectednumberofhits,level,col6,pvals_adjusted[a],'MF',sep=";")
	    cat("\n")
    }
    rm(targeted_goterm)
    a=a+1
  }
  ,file=paste(outputPrefix,"/hyperesults_mf.txt",sep=""))

#################### CELLULAR COMPOPENT #################

Go_univ=univers_cc
Go_ech=echantillon_cc
Go_uniq=Go_univ[!duplicated(Go_univ),]

#####PHYPER TEST WITH BONFERRONI CORRECTION#####
pvals=""
for(i in 1:length(Go_uniq))
{
  targeted_goterm=toString(Go_uniq[i]) #tostring in order to request value
  pvals[i]<-phyper(length(which(Go_ech==targeted_goterm)),length(which(Go_univ==targeted_goterm)),length(Go_univ[,1])-length(which(Go_univ==targeted_goterm)),length(Go_ech[,1]),lower.tail = FALSE) 
}
pvals2<-as.vector(pvals)
pvals_adjusted<-p.adjust(pvals2, method="bonferroni", n=length(pvals))
a=1
######PHYPER OK, NOW ENHANCE PHYPER RESULTS#####
capture.output(
  for(i in 1:length(Go_uniq))
  {
    targeted_goterm=toString(Go_uniq[i]) #tostring in order to request value
    expectednumberofhits=(((length(which(Go_univ==targeted_goterm)))/(length(Go_univ[,1])))*length(Go_ech[,1]))
    trm=Term(targeted_goterm)
    ##########################################################################
    #RETRIEVE GO-TERM LEVEL : by parents "walking" method (is-a , part-of)
    ##########################################################################
    #First request
    parents=get(targeted_goterm,GOCCPARENTS) 
    #If parents founded-->next DAG level exists.So repeat request with the parents of the next DAG level
    if(parents[1]!=parent)
    {
      level=1 #level in the DAG
      tryCatch({
        for(i in 1:1000) #Arbitrary.
        {
          parents=get(parents[1],GOCCPARENTS) #[1] because prioritize for "is-a"
          level=level+1
          if(parents[1]==parent) stop()#If node ="all" , all go three processed,no more DAG levels.Stop the requests loop
        }
      }, error=function(e){cat()})#Nothing in message error in order to dont print msgerros in results file
    }
    #If no parents founded-->term is on highest level. This is Aspect term-->lvl
    else
    {
      level=1
    }
    #Print results foreach go-term
    if(length(which(Go_ech==targeted_goterm))!=0) #if go-term is present in sample
    {
	    col3=as.character(length(which(Go_ech==targeted_goterm)))
	    col6=as.character(phyper(length(which(Go_ech==targeted_goterm)),length(which(Go_univ==targeted_goterm)),length(Go_univ[,1])-length(which(Go_univ==targeted_goterm)),length(Go_ech[,1]),lower.tail = FALSE))
	    cat(targeted_goterm,trm,col3,expectednumberofhits,level,col6,pvals_adjusted[a],'CC',sep=";")
	    cat("\n")
    }
    rm(targeted_goterm)
    a=a+1
  }
  ,file=paste(outputPrefix,"/hyperesults_cc.txt",sep=""))

#################### BIOLOGICAL PROCESS #################

Go_univ=univers_bp
Go_ech=echantillon_bp
Go_uniq=Go_univ[!duplicated(Go_univ),]

#####PHYPER TEST WITH BONFERRONI CORRECTION#####
pvals=""
for(i in 1:length(Go_uniq))
{
  targeted_goterm=toString(Go_uniq[i]) #tostring in order to request value
  pvals[i]<-phyper(length(which(Go_ech==targeted_goterm)),length(which(Go_univ==targeted_goterm)),length(Go_univ[,1])-length(which(Go_univ==targeted_goterm)),length(Go_ech[,1]),lower.tail = FALSE) 
}
pvals2<-as.vector(pvals)
pvals_adjusted<-p.adjust(pvals2, method="bonferroni", n=length(pvals))
a=1

######PHYPER OK, NOW ENHANCE PHYPER RESULTS#####
capture.output(
  for(i in 1:length(Go_uniq))
  {
    targeted_goterm=toString(Go_uniq[i]) #tostring in order to request value
    expectednumberofhits=(((length(which(Go_univ==targeted_goterm)))/(length(Go_univ[,1])))*length(Go_ech[,1]))
    trm=Term(targeted_goterm)
    ##########################################################################
    #RETRIEVE GO-TERM LEVEL : by parents "walking" method (is-a , part-of)
    ##########################################################################
    #First request
    parents=get(targeted_goterm,GOBPPARENTS) 
    #If parents founded-->next DAG level exists.So repeat request with the parents of the next DAG level
    if(parents[1]!=parent)
    {
      level=1 #level in the DAG
      tryCatch({
        for(i in 1:1000) #Arbitrary.
        {
          parents=get(parents[1],GOBPPARENTS) #[1] because prioritize for "is-a"
          level=level+1
          if(parents[1]==parent) stop()#If node ="all" , all go three processed,no more DAG levels.Stop the requests loop
        }
      }, error=function(e){cat()})#Nothing in message error in order to dont print msgerros in results file
    }
    #If no parents founded-->term is on highest level. This is Aspect term-->lvl
    else
    {
      level=1
    }
    #Print results foreach go-term
    if(length(which(Go_ech==targeted_goterm))!=0) #if go-term is present in sample
    {
	    col3=as.character(length(which(Go_ech==targeted_goterm)))
	    col6=as.character(phyper(length(which(Go_ech==targeted_goterm)),length(which(Go_univ==targeted_goterm)),length(Go_univ[,1])-length(which(Go_univ==targeted_goterm)),length(Go_ech[,1]),lower.tail = FALSE))
	    cat(targeted_goterm,trm,col3,expectednumberofhits,level,col6,pvals_adjusted[a],'BP',sep=";")
	    cat("\n")
	}
    rm(targeted_goterm)
    a=a+1
  }
  ,file=paste(outputPrefix,"/hyperesults_bp.txt",sep=""))
