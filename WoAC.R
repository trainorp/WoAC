############ Prereqs ############
args<-commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=FALSE)
library(methods)
library(nnet)
library(cvTools)
library(doParallel)
library(dplyr)

############ Bootstrap ############
if(as.logical(as.character(args[2]))){
  set.seed(as.numeric(args[1]))
  boot<-sample(1:nrow(DF),replace=TRUE)
  set.seed(as.numeric(args[1])+1000)
  randSub<-sample(1:ncol(DF),size=floor(ncol(DF)/4))
  DF<-DF[boot,randSub]
  group<-group[boot,]
}

############ Population class ############
cl<-makeCluster(2)
registerDoParallel(cl)

setClass("population",slots=list(varInclude="list",age="integer",cost="numeric"))
cost<-function(x){
  return(unlist(lapply(x,FUN=function(y){
    if(length(table(y))<2){
      return(Inf)
    }else{
      DF2<-cbind(group=group$group,as.data.frame(DF[,y]))
      DF2$group<-factor(DF2$group)
      cv1<-cvFolds(nrow(DF2),K=10,R=10)
      mis<-foreach(j=1:10,.combine=c,.packages=c("nnet"),.inorder=FALSE) %dopar% {
        mis2<-c()
        cv2<-data.frame(samp=cv1$subsets[,j],fold=cv1$which)
        for(i in 1:10){
          train<-cv2$samp[cv2$fold!=i]
          test<-cv2$samp[cv2$fold==i]
          
          pred1<-log(predict(multinom(group~.,data=DF2[train,],
                                      MaxNWts=10000,maxit=40,reltol=.0001,trace=FALSE),
                             newdata=DF2[test,],type="probs"))
          pred1[is.infinite(pred1)]<-log(1e-200)
          mm1<-model.matrix(~-1+DF2[test,]$group)
          mis2<-c(mis2,-1/nrow(mm1)*sum(mm1*pred1))
        }
        mean(mis2)
      }
      return(mean(mis))
    }
  })))
}

setGeneric("minCost",function(x,n=1) standardGeneric("minCost"))
setMethod("minCost","population",
          function(x,n){
            return(order(x@cost)[1:n])
          }
)

setGeneric("relFitness",function(x,inverse=FALSE) standardGeneric("relFitness"))
setMethod("relFitness","population",
          function(x,inverse){
            fits<-x@cost
            ecdf1<-ecdf(fits)(fits)
            if(!inverse){
              return(1-ecdf1+min(ecdf1))
            }else{
              return(ecdf1)
            }
          }
)

setGeneric("recombination",function(x,beta=3) standardGeneric("recombination"))
setMethod("recombination","population",
          function(x,beta){
            relFit<-relFitness(x)
            tinderProb<-qbeta(relFit,shape1=1,shape2=beta)
            tinder<-which(sapply(tinderProb,FUN=function(y) as.logical(rbinom(n=1,size=1,prob=y))))
            weds<-matrix(NA,nrow=0,ncol=2)
            while(length(tinder)>1){
              pair<-sample(tinder,size=2)
              weds<-rbind(weds,pair)
              tinder<-tinder[!tinder %in% pair]
            }
            if(nrow(weds)>1){
              children<-list()
              for(i in 1:nrow(weds)){
                mom<-x@varInclude[[weds[i,][1]]]
                dad<-x@varInclude[[weds[i,][2]]]
                # Two random crossover ends:
                cross<-sample(1:(length(mom)-1),2)
                cross<-cross[order(cross)]
                children<-c(children,
                            list(c(dad[1:cross[1]],mom[(cross[1]+1):cross[2]],dad[(cross[2]+1):length(dad)])))
              }
              x@varInclude<-c(x@varInclude,children)
              x@age<-c(x@age,integer(length(children)))
              x@cost<-c(x@cost,cost(children))
            }
            return(x)
          }
)

setGeneric("mutation",function(x) standardGeneric("mutation"))
setMethod("mutation","population",
          function(x){
            relFit<-relFitness(x,inverse=TRUE)
            invFit<-(relFit-min(relFit))/10
            for(j in 1:length(invFit)){
              phi<-invFit[j]
              piT<-5/ncol(DF)
              piF<-1-piT
              aTT<-piT+piF*(1-phi)
              aTF<-1-aTT
              aFT<-phi-aTF
              aFF<-1-aFT
              A<-matrix(c(aTT,aTF,aFT,aFF),byrow=TRUE,nrow=2)
              a<-t(A)
              rownames(a)<-colnames(a)<-c("TRUE","FALSE")
              
              trues<-which(x@varInclude[[j]])
              TSwitch<-sample(trues,size=min(rbinom(1,size=length(trues),a['FALSE','TRUE']),
                                             length(trues)))
              x@varInclude[[j]][TSwitch]<-FALSE
              
              falses<-which(!x@varInclude[[j]])
              FSwitch<-sample(falses,size=min(rbinom(1,size=length(falses),a['TRUE','FALSE']),
                                              length(falses)))
              x@varInclude[[j]][FSwitch]<-TRUE
            }
            x@cost<-cost(x@varInclude)
            return(x)
          }
)

setGeneric("death",function(x,popSize) standardGeneric("death"))
setMethod("death","population",
          function(x,popSize){
            costs<-x@cost
            infs<-which(is.infinite(costs))
            if(length(infs)>0){
              x@varInclude<-x@varInclude[-infs]
              x@age<-x@age[-infs]
              x@cost<-x@cost[-infs]
            }
            
            relFit<-relFitness(x)
            ageFit<-ecdf(x@age)(x@age)
            totalFit<-3*relFit+ageFit
            killN<-length(x@varInclude)-popSize
            if(killN>0){
              x@varInclude<-x@varInclude[-order(totalFit)[1:killN]]
              x@age<-x@age[-order(totalFit)[1:killN]]
              x@cost<-x@cost[-order(totalFit)[1:killN]]
            }
            return(x)
          }
)

setGeneric("migration",function(x,varProp,migrationP=0.1) standardGeneric("migration"))
setMethod("migration","population",
          function(x,varProp,migrationP){
            migrants<-ceiling(length(x@varInclude)*migrationP)
            migrants<-birth(migrants,varProp)
            x@varInclude<-c(x@varInclude,migrants@varInclude)
            x@age<-c(x@age,migrants@age)
            x@cost<-c(x@cost,migrants@cost)
            return(x)
          }
)

birth<-function(births,varProp)
{
  success<-rbinom(births,size=ncol(DF),prob=varProp)
  varInclude<-lapply(success,FUN=function(x){
    y<-rep(FALSE,ncol(DF))
    y[sample(1:length(y),x)]<-TRUE
    return(y)
  })
  return(new("population",varInclude=varInclude,age=integer(length(varInclude)),
             cost=cost(varInclude)))
}

generation<-function(popSize=250,varProp=5/ncol(DF),migrationP=.1,iterations=200)
{
  progress<-c()
  births<-popSize
  pop0<-birth(births,varProp)
  progress<-rbind(progress,data.frame(action="start",size=length(pop0@varInclude),which=minCost(pop0),val=min(pop0@cost)))
  for(j in 1:iterations){
    pop0@age<-pop0@age+as.integer(1)
    pop0<-recombination(pop0)
    progress<-rbind(progress,data.frame(action="recombination",size=length(pop0@varInclude),which=minCost(pop0),val=min(pop0@cost)))
    pop0<-mutation(pop0)
    progress<-rbind(progress,data.frame(action="mutation",size=length(pop0@varInclude),which=minCost(pop0),val=min(pop0@cost)))
    pop0<-migration(pop0,varProp,migrationP=migrationP)
    pop0<-death(pop0,popSize)
    progress<-rbind(progress,data.frame(action="migration",size=length(pop0@varInclude),which=minCost(pop0),val=min(pop0@cost)))
    cat("j is:",j,"\n")
  }
  list(pop=pop0,vars=names(DF),trace=progress)
}

########## Generate Crowds ##########
out<-list()
ptm<-proc.time()
for(i in 1:10){
  cat("i is:",i,"\n")
  out[[i]]<-generation()
}
print(proc.time()-ptm)
stopCluster(cl)

assign(paste0("out",args[1]),out)
save(list=paste0("out",args[1]),file=paste0("out",args[1],".RData"))