#QBL_V5 is an additional version of QBLstrat (https://github.com/mxw010/LBL). QBL_V5 is built under the same framework as QBLstrat, but implements normal priors for PC effects, instead of Laplace priors. The prior for PC effects is i.i.d N(0, sigma_a^2), and sigma_a^2 is assigned with an Inverse-Gamma(randomA_sigmas_a, randomA_sigmas_a) hyper-prior.

#Below are some new arguments for QBL_V5:

#' @param R The PC data to adjust for population stratification. The data set should consist of n rows and n.pc columns, where n.pc is the number of PCs. Column names should be specified. PC effects are assigned with normnal priors.
#' @param start.gamma Starting value of all regression coefficients for PCs, default is 0.01.
#' @param randomA_sigmas Starting value of the \eqn{\sigma_a^2} parameter; default is 0.7.
#' @param randomA_sigmas_a First hyperparameter of the Inverse-Gamma(\eqn{randomA_sigmas_a}, \eqn{randomA_sigmas_b}) prior for the variance of the normal PC effect, \eqn{\sigma_a^2}. The default value is 1.
#' @param randomA_sigmas_b Second hyperparameter of Inverse-Gamma prior described above; default is 1.

#' Return a list with the following components:
#' \describe{
#' \item{BF}{Bayes Factors for haplotypes and covariates (e.g., age and sex).}
#' \item{beta}{The coefficient estimates for haplotypes and covariates.}
#' \item{gamma}{The coefficient estimates for PCs.}
#' \item{CI.beta}{The 95% credible intervals for haplotypes and covariates.}
#' \item{CI.gamma}{The 95% credible intervals for PCs.}
#' \item{CI.lambda}{The 95% credible intervals for \eqn{\lambda}.}
#' \item{CI.sigmas}{The 95% credible intervals for \eqn{\sigma^2}.}
#' \item{CI.sigmasA}{The 95% credible intervals for \eqn{\sigma_a^2}.}
#' \item{CI.D}{The 95% credible intervals for D.}
#' \item{freq}{The estimated haplotype frequencies.}
#'}
#'



QBL_V5<-function(dat, numSNPs=5, allelic=TRUE, baseline = "missing", interaction=F,cov=1, pooling.tol=0,zero.tol="missing",
                      a = 20, b = 20, start.beta = 0.01, start.gamma = 0.01, lambda = 1, D = 0, seed = NULL, e = 0.1, burn.in = 20000,
                      num.it = 50000,sigmas=1,asig=2,bsig=1,SDD,beta.mu.0=0,beta.sigma.0=1,dirres.resl,l,setting, CC=1000, 
								      randomA_sigmas=0.7,randomA_sigmas_a=1,randomA_sigmas_b=1,R, dumpConvergence=F,n.chains=1,plot.interval=NULL,
								      checkGamma=c(1:5),checkCores=3)

{
	require(coda)
	require(smoothmest)
	
  R_vec<-c(t(R))
  ### Prior settings
  cat("\n Prior settings before qLBL C program\n")
  cat("a is", a, "\n")
  cat("b is", b, "\n")
  cat("asigma is", asig, "\n")
  cat("bsigma is", bsig, "\n")
  cat("beta.mu.0 is", beta.mu.0, "\n")
  cat("bsigma is", beta.sigma.0, "\n")
  cat("\n End of Prior settings before qLBL C program\n")
  #########
  
  if(zero.tol=="missing"){
        haplos.new.list<-pre.hapassoc(dat, numSNPs=numSNPs, pooling.tol=pooling.tol, allelic=allelic,verbose=F)
  }else{
        haplos.new.list<-pre.hapassoc(dat, numSNPs=numSNPs, pooling.tol=pooling.tol,zero.tol=zero.tol, allelic=allelic,verbose=F)
  }
  
  cat("pre.hapassoc")
  
  
  #  cat("\n","initial frequence names",  names(haplos.new.list$initFreq), "\n")
  #  cat("\n","initial frequence",  haplos.new.list$initFreq, "\n")  
  # freq.initial.test[l,]<-names(haplos.new.list$initFreq)
  
  
  if (is.null(seed)==FALSE) set.seed(seed)
  haplos.names <- names(haplos.new.list$initFreq)
  freq <- haplos.new.list$initFreq
  if (baseline=="missing")
  {
    baseline <- haplos.names[which.max(freq)]
  }
  
  
  # cat("\n","maximum frequency from pre.hapassoc",  haplos.names[which.max(freq)], "\n")
  # freq.initial.max[l]<-haplos.names[which.max(freq)]
  # # 
  # if(freq.initial.max[l]==baseline)
  # {
  #   mmm=mmm+1
  # }
  
  
  column.subset <- colnames(haplos.new.list$haploDM) != baseline
  n.cov<-dim(haplos.new.list$nonHaploDM)[2] - 1
  
  
  #cov.data<-haplos.new.list$nonHaploDM[,-1]
  
  #hdat <- cbind(haplos.new.list$nonHaploDM[,1], haplos.new.list$haploDM[, column.subset], cov.data)
  #colnames(hdat) <- c(colnames(haplos.new.list$nonHaploDM)[1], colnames(haplos.new.list$haploDM[, column.subset]), 
  #                    names(haplos.new.list$nonHaploDM)[-1])
  
  
  
  
  
  if (cov!=0){
    cov.data<-haplos.new.list$nonHaploDM[,-1]     ## it is 2 previous ###
    names.cov<-names(haplos.new.list$nonHaploDM)[-1]  ## it is 2 previous ###
    
    if(interaction==T)
    {
      if(cov>1){
        
        cov.data.int <- cov.data ## create and add dummy interaction columns
        cov.data.int<-cbind(haplos.new.list$haploDM[, column.subset]*cov.data[,1],haplos.new.list$haploDM[, column.subset]*cov.data[,2],cov.data.int)
        
        
        
        t2<-dim(cbind(haplos.new.list$haploDM[, column.subset]*cov.data[,1],haplos.new.list$haploDM[, column.subset]*cov.data[,2]))[2]
        colnames(cov.data.int)[1:t2]<-c(paste(names(haplos.new.list$haploDM[, column.subset,drop=F]), names.cov[1], sep=""),paste(names(haplos.new.list$haploDM[, column.subset,drop=F]), names.cov[2], sep=""))
        colnames(cov.data.int)[c(t2+1,t2+2)]<-names.cov
        
        hdat <- cbind(haplos.new.list$nonHaploDM[,1], haplos.new.list$haploDM[, column.subset], cov.data.int)     
        colnames(hdat) <- c(colnames(haplos.new.list$nonHaploDM)[1], colnames(haplos.new.list$haploDM[, column.subset,drop=F]), colnames(cov.data.int))
        
        
        
        
      }else{ #cov=1
        
        cov.data.int <- cov.data ## create and add dummy interaction columns
        
        cov.data.int<-cbind(haplos.new.list$haploDM[, column.subset]*cov.data,cov.data.int)
        
        t2<-dim(haplos.new.list$haploDM[, column.subset,drop=F]*cov.data)[2]
        
        #t1<-dim(dum.data.int)[2]
        #t<-t1-t2+1
        colnames(cov.data.int)[1:t2]<-paste(names(haplos.new.list$haploDM[, column.subset,drop=F]), names.cov, sep="")
        colnames(cov.data.int)[t2+1]<-names.cov
        
        
        hdat <- cbind(haplos.new.list$nonHaploDM[,1], haplos.new.list$haploDM[, column.subset], cov.data.int)     
        colnames(hdat) <- c(colnames(haplos.new.list$nonHaploDM)[1], colnames(haplos.new.list$haploDM[, column.subset,drop=F]), colnames(cov.data.int))
      }
      
      
      
    }else{ #interaction=F
      
      
      
      hdat <- cbind(haplos.new.list$nonHaploDM[,1], haplos.new.list$haploDM[, column.subset], cov.data)
      #hdat <- cbind(haplos.new.list$nonHaploDM[,1], haplos.new.list$haploDM[, column.subset])
      colnames(hdat) <- c(colnames(haplos.new.list$nonHaploDM)[1], colnames(haplos.new.list$haploDM[, column.subset,drop=F]), 
                          names(haplos.new.list$nonHaploDM)[-1])
      #colnames(hdat) <- c(colnames(haplos.new.list$nonHaploDM)[1], colnames(haplos.new.list$haploDM[, column.subset]))
    }
    
  }else{ #cov=0
    cat("no cov")
    hdat <- cbind(haplos.new.list$nonHaploDM[,1],haplos.new.list$haploDM[, column.subset])
    colnames(hdat) <- c(colnames(haplos.new.list$nonHaploDM)[1],colnames(haplos.new.list$haploDM[, column.subset,drop=F]))
    
  }
  
  hap.name=colnames(hdat[,2:dim(haplos.new.list$haploDM)[2]])
  
  
  #### Check what hdat is
  cat("head of hdat is\n")
  #print(head(hdat))
  #cat("class of hdat is ",class(hdat),"\n")
  #stop("checking point")
  ####
  
  
  
  #############
  #n.cov.dum<-dim.dum
  ID <- haplos.new.list$ID
  N <- sum(haplos.new.list$wt)
  
  if(is.null(ncol(R))){
      q <- 1
   }else{
      q <- ncol(R)
   }

  y <- as.numeric(hdat[, 1]) #trait value
  x <- data.matrix(hdat[, -1,drop=F])  #hap data, include baseline
  colnames(x)<-NULL
  
  
  freq.new<-freq[names(haplos.new.list$haploDM[, column.subset,drop=F])]
  freq.new[length(freq.new)+1]<-freq[baseline]
  freq.new<-as.vector(freq.new)
  num.haplo.id<-as.vector(table(ID))
  x.length<-as.integer(dim(x)[2]) 
  
  h.length<-dim(haplos.new.list$haploDM)[2]
  freq.data<-haplos.new.list$haploDM[, column.subset,drop=F] 
  freq.num<-matrix(rep(NA,2*(dim(freq.data)[1])),ncol=2)
  for(i in 1:dim(freq.data)[1])
  {
    for(j in 1:dim(freq.data)[2])
    {
      if(freq.data[i,j]==2)
      {
        freq.num[i,1]=j
        freq.num[i,2]=j
        break
      }
      if(freq.data[i,j]==1)
      {
        if(is.na(freq.num[i,1])==T) freq.num[i,1]=j
        else
        {
          freq.num[i,2]=j
          break
        }
      }
    }
    if(is.na(freq.num[i,1])==T)
    {
      freq.num[i,1]=dim(freq.data)[2]+1
      freq.num[i,2]=dim(freq.data)[2]+1
    }
    if(is.na(freq.num[i,2])==T) freq.num[i,2]=dim(freq.data)[2]+1
  }
  
  
  
  
  
  
  ########TEST######
  #name.cs=names(table(haplos.new.list$nonHaploDM[,2]))
  #name.cs=name.cs[name.cs!=baseline.cov]
  #len.E=length(name.cs)+1
  #Num.E=rep(NA,len.E)
  #Num.E[1]=0
  #for(i in 2:len.E)
  #{
  #  Num.E[i]=length(ddat[ddat[,2]==name.cs[i-1],2])
  #}
  #Num.E[1]=N-sum(Num.E)
  ########
  
  beta=rep(start.beta, x.length+1) #initial all beta_j to a small number
  gamma=rep(start.gamma, q)
  beta.out<-numeric((num.it-burn.in)*(x.length+1)) # num of sampled beta's
  gamma.out<-numeric((num.it-burn.in)*q)
  lambda.out<-numeric(num.it-burn.in) # num of sampled lambda's(prior of beta)
  sigmas.out<-numeric(num.it-burn.in)
  randomA_sigmas.out<-numeric(num.it-burn.in)
  freq.out<-numeric((num.it-burn.in)*(h.length)) #num of sampled f's
  D.out<-numeric(num.it-burn.in) # number of 'd'
  #beta=rep(start.beta, x.length+1) #initial all beta_j to a small number
  #beta.out<-numeric((num.it-burn.in)*(x.length+1)) # num of sampled beta's
  
  up.xz=numeric(x.length*N)
  
  BDIC.out<-numeric(num.it-burn.in) # BDIC
  
  
  ##### Before put into the program, check those initial values
  
  cat("\n Before put into the program, check those initial values is\n")
  cat("Head of x is\n")
  print(head(x))
  cat("Dimension of x is \n")
  print(dim(x)[2])
  #cat("Number of people in design matrix x is ", dim(x)[1],"\n")
  #cat("first several responses\n")
  #print(head(y))
  #cat("Total number of people in original design matrix is",N,"\n")
  #cat("Head of what haplotype pairs each person has\n")
  #head(num.haplo.id)
  #cat("x.length is equal to",x.length,"\n")
  #cat("total number of possible haplotypes is",h.length,"\n")
  #cat("Head of freq.num is\n")
  #print(head(freq.num))
  #cat("freq.new is\n")
  #print(freq.new)        
  #cat("The initial value of D is",D,"\n")        
  #cat("The initial value of sigmas is ",sigmas,"\n")
  #cat("The initial value of beta is \n")
  #print(beta)        
  #cat("The initial value of a is",a,"\n")
  #cat("The initial value of b is",b,"\n")
  #cat("The initial value of asig is",asig,"\n")
  #cat("The initial value of bsig is",bsig,"\n")
  #cat("Value of num.it is",num.it,"\n")         
  #cat("Value of burn.in is",burn.in,"\n")         
  #cat("number of number of covariate is",n.cov,"\n")        
  
  #cat("Initial value of lambda is",lambda,"\n")
  #cat("Initial values beta.out is\n");
  #print(beta.out)
  #cat("First lambda.out is ", lambda.out[0],"\n");
  #cat("First sigmas.out is ", sigmas.out[0],"\n");
  #cat("Initial freq.out is \n");
  #print(freq.out)
  #cat("Initial D.out is ", D.out,"\n\n");
  
  cat("Prior mean of beta.0 is ",beta.mu.0,"\n")
  cat("Prior variance of beta.0 is ",beta.sigma.0,"\n")
  
  
  
  
  
  cat("End before put into the program, check those initial values is\n")
  #####
  
  
  dyn.load("QBL_V5.so")
  #dyn.load("LBLq_han_d_3.so")
  #dyn.load("LBLq_han_d_0.so")
  
  #start from here
  #ZL: check convergence
  if(dumpConvergence==T){
    
    #XF: plot.interval
    if(is.null(plot.interval)==T){
      plot.interval<-(num.it-burn.in)%/%300
    }
    
    my.draws<-vector("list",n.chains)
    beta.matrix<-vector("list",n.chains)
    gamma.matrix<-vector("list",n.chains)
    lambda.vector<-matrix(numeric(n.chains*(num.it-burn.in)),ncol=n.chains)
    sigmas.vector<-matrix(numeric(n.chains*(num.it-burn.in)),ncol=n.chains)
    randomA_sigmas.vector<-matrix(numeric(n.chains*(num.it-burn.in)),ncol=n.chains)
    
    
    
    name1<-colnames(hdat[, -1,drop=F])
    if(q==1){
        names.gamma<-"PC1" #if R n by 1, name is PC1
    }else{
        names.gamma<-colnames(R)
    }
    names<-c("Intercept",name1,names.gamma,"lambda","sigmas","randomA_sigmas")          #XF: changes made
    

    multiRun<-function(c){
    	
    	start.lambda<-rgamma(1,a,b)
    	start.beta<-rdoublex(1,mu=0,lambda=1/start.lambda)
    	start.randomA_sigmas<-rgamma(1,shape=randomA_sigmas_a,rate=randomA_sigmas_b)
    	start.randomA_sigmas<-1/start.randomA_sigmas
    	start.gamma<-rnorm(1,0,sqrt(start.randomA_sigmas))
    	beta<-rep(start.beta, x.length+1)
    	gamma<-rep(start.gamma, q)
    	
    	out<-.C("GXE_mcmc", x=as.double(x),
    					n=as.integer(dim(x)[1]), as.double(y),
    					as.integer(N),as.integer(q), as.integer(num.haplo.id),
    					as.integer(x.length), as.integer(h.length),
    					as.integer(freq.num), as.double(freq.new),
    					as.double(D), as.double(sigmas),as.double(beta), as.double(gamma), as.double(a), as.double(b),
    					as.double(asig),as.double(bsig),  as.double(start.lambda), as.integer(num.it),
    					as.integer(burn.in), beta.out=as.double(beta.out), gamma.out=as.double(gamma.out),
    					lambda.out=as.double(lambda.out),  sigmas.out=as.double(sigmas.out),
    					randomA_sigmas.out=as.double(randomA_sigmas.out),  freq.out=as.double(freq.out),D.out=as.double(D.out),as.integer(n.cov),up.xz=as.double(up.xz),as.integer(CC),BDIC.out=as.double(BDIC.out),
    					as.double(start.randomA_sigmas),as.double(randomA_sigmas_a),as.double(randomA_sigmas_b),as.double(R_vec))
    	

    	beta.outn<-matrix(out$beta.out,nrow=num.it-burn.in, byrow=TRUE)
    	gamma.outn<-matrix(out$gamma.out,nrow=num.it-burn.in, byrow=TRUE)
    	lambda.outn<-matrix(out$lambda.out,nrow=num.it-burn.in, byrow=TRUE)
    	sigmas.outn<-matrix(out$sigmas.out,nrow=num.it-burn.in, byrow=TRUE)
    	randomA_sigmas.outn<-matrix(out$randomA_sigmas.out,nrow=num.it-burn.in, byrow=TRUE)
    	
    	beta_lambda_result<-cbind(beta.outn,gamma.outn,lambda.outn,sigmas.outn,randomA_sigmas.outn)
    	colnames(beta_lambda_result)<-names
    	
    	return(beta_lambda_result)
    	
    }
    
    #parallel run
    #registerDoParallel(checkCores)
    #checkConverg<-foreach(c=1:n.chains) %dopar% {
    #	multiRun(c)
    #}
    #stopImplicitCluster()
    
    #lapply
    checkConverg<-lapply(c(1:n.chains),FUN=multiRun)
    
    t<-ncol(checkConverg[[1]])
    for(i in 1:n.chains){
    	
    	my.draws[[i]]<-mcmc(checkConverg[[i]])
    	beta.matrix[[i]]<-checkConverg[[i]][,1:(t-3-q)]
    	gamma.matrix[[i]]<-checkConverg[[i]][,(t-2-q):(t-3)]
    	lambda.vector[,i]<-checkConverg[[i]][,t-2]
    	sigmas.vector[,i]<-checkConverg[[i]][,t-1]
    	randomA_sigmas.vector[,i]<-checkConverg[[i]][,t]
    }

    mh.list <- mcmc.list(my.draws)
    
    if(n.chains==1){
      D_Raf<-raftery.diag(my.draws[[1]],q=0.025,r=0.005,s=0.95)
      num.beta<-dim(beta.matrix[[1]])[2]
      #pos.lambda<-num.beta+1       #dim(beta.matrix[[1]])[2] is num.beta
      acceptance.rate<-1-rejectionRate(my.draws[[1]])
      
      #Convergence Numeric Summary
      sink("QBLDiagnosticNumericalSummary.txt")
      cat("Gelman Diagnostic needs at least 2 chains!\n\n\n")
      cat("Raftery Diagnostic \n")
      print(D_Raf)
      cat("Acceptance Rate \n")
      cat(acceptance.rate,"\n")
      sink()
      
      #Convergence Graphical Summary
      pdf(file = "QBLDiagnosticPlot.pdf")
      
      par(mfrow=c(3,2))
      
      #draw trace plot for beta
      seq.get<-seq(1,length(lambda.vector[,1]),plot.interval)
      
      for (i in 1: num.beta){
        plot(seq.get,beta.matrix[[1]][seq.get,i],type="l",col=1,ylim=c(min(beta.matrix[[1]][seq.get,i])-0.01,max(beta.matrix[[1]][seq.get,i])+0.01),
              ylab=paste("beta",i,sep=""),xlab="iteration")
        #plot(seq.get,beta.matrix[[1]][seq.get,i],type="l",col=1,ylim=c(-2,2),
        #     ylab=paste("beta",i,sep=""),xlab="iteration")
        abline(h=0,col=(n.chains+1))
      }
      
      #draw trace plot for gamma
      for (i in 1: q){
        plot(seq.get,gamma.matrix[[1]][seq.get,i],type="l",col=1,ylim=c(min(gamma.matrix[[1]][seq.get,i])-0.01,max(gamma.matrix[[1]][seq.get,i])+0.01),
             ylab=paste("gamma",i,sep=""),xlab="iteration")
        #plot(seq.get,gamma.matrix[[1]][seq.get,i],type="l",col=1,ylim=c(-2,2),
        #     ylab=paste("gamma",i,sep=""),xlab="iteration")
        abline(h=0,col=(n.chains+1))
      }
      
      #draw trace plot for lambda
      plot(seq.get,lambda.vector[seq.get,1],type="l",col=1,ylim=c(min(lambda.vector[seq.get,1])-0.01,max(lambda.vector[seq.get,1])+0.01),
           ylab="lambda",xlab="iteration")
      abline(h=0,col=(n.chains+1))
      
      #draw trace plot for sigmas
      plot(seq.get,sigmas.vector[seq.get,1],type="l",col=1,ylim=c(min(sigmas.vector[seq.get,1])-0.01,max(sigmas.vector[seq.get,1])+0.01),
               ylab="sigmas",xlab="iteration")
      abline(h=0,col=(n.chains+1))
      
      #draw trace plot for randomA_sigmas
      plot(seq.get,randomA_sigmas.vector[seq.get,1],type="l",col=1,ylim=c(min(randomA_sigmas.vector[seq.get,1])-0.01,max(randomA_sigmas.vector[seq.get,1])+0.01),
               ylab="randomA_sigmas",xlab="iteration")
      abline(h=0,col=(n.chains+1))
      dev.off()
      
    }else if(n.chains>1){
      D_Gel<-gelman.diag(mh.list)
      
      D_Raf<-vector("list",n.chains)
      acceptance.rate<-vector("list",n.chains)
      
      num.beta<-dim(beta.matrix[[1]])[2]
      #pos.lambda<-num.beta+1 #position of lambda
      
      for (i in 1:n.chains){
        D_Raf[[i]]<-raftery.diag(my.draws[[i]],q=0.025,r=0.005,s=0.95)
        acceptance.rate[[i]]<-1-rejectionRate(my.draws[[i]])   #acceptance rate of one MCMC
      }
      
      #Convergence Numeric Summary
      sink("QBLDiagnosticNumericalSummary.txt")
      cat("Gelman Diagnostic \n")
      print(D_Gel)
      cat("\n\n\n")
      
      for (i in 1:n.chains){
        cat("Chain ",i,"'s Result \n",sep="")
        cat("Raftery Diagnostic \n")
        print(D_Raf[[i]])
        cat("Acceptance Rate \n")
        print(acceptance.rate[[i]])
        cat("\n\n\n")
      }
      sink()
      
      #Convergence Graphical Summary
      pdf(file = "QBLDiagnosticPlot.pdf")
      
      # draw Gelman plot
      gelman.plot(mh.list)
      
      par(mfrow=c(3,2))
      
      #draw trace plot for beta
      seq.get<-seq(1,length(lambda.vector[,1]),plot.interval)
      
      for (i in 1: num.beta){
        plot(seq.get,beta.matrix[[1]][seq.get,i],type="l",col=1,ylim=c(min(beta.matrix[[1]][seq.get,i])-0.5,max(beta.matrix[[1]][seq.get,i])+0.5),
             ylab=names[i],xlab="iteration")
        # plot(seq.get,beta.matrix[[1]][seq.get,i],type="l",col=1,ylim=c(-2,2),
        #      ylab=paste("beta",i,sep=""),xlab="iteration")
        for (j in 2:n.chains){
          lines(seq.get,beta.matrix[[j]][seq.get,i],col=j)
        }
        abline(h=0,col=(n.chains+1))
      }
      
      #draw trace plot for gamma: checkGamma
      for (i in checkGamma){
      	plot(seq.get,gamma.matrix[[1]][seq.get,i],type="l",col=1,ylim=c(min(gamma.matrix[[1]][seq.get,i])-0.5,max(gamma.matrix[[1]][seq.get,i])+0.5),
        ylab=names.gamma[i],xlab="iteration")
        # plot(seq.get,beta.matrix[[1]][seq.get,i],type="l",col=1,ylim=c(-2,2),
        #      ylab=paste("beta",i,sep=""),xlab="iteration")
        for (j in 2:n.chains){
          lines(seq.get,gamma.matrix[[j]][seq.get,i],col=j)
        }
        abline(h=0,col=(n.chains+1))
      }
      
      #draw trace plot for lambda
      plot(seq.get,lambda.vector[seq.get,1],type="l",col=1,ylim=c(min(lambda.vector[seq.get,1])-0.5,max(lambda.vector[seq.get,1])+0.5),ylab="lambda"
           ,xlab="iteration")
      for (i in 2: n.chains){
        lines(seq.get,lambda.vector[seq.get,i],col=i,ylab="lambda")
      }
      abline(h=0,col=(n.chains+1))
      
      #draw trace plot for sigmas
      plot(seq.get,sigmas.vector[seq.get,1],type="l",col=1,ylim=c(min(sigmas.vector[seq.get,1])-0.5,max(sigmas.vector[seq.get,1])+0.5),
      		 ylab="sigmas",xlab="iteration")
      for (i in 2: n.chains){
      	lines(seq.get,sigmas.vector[seq.get,i],col=i,ylab="sigmas")
      }
      abline(h=0,col=(n.chains+1))
      
      #draw trace plot for randomA_sigmas
      plot(seq.get,randomA_sigmas.vector[seq.get,1],type="l",col=1,ylim=c(min(randomA_sigmas.vector[seq.get,1])-0.5,max(randomA_sigmas.vector[seq.get,1])+0.5),
      		 ylab="randomA_sigmas",xlab="iteration")
      for (i in 2: n.chains){
      	lines(seq.get,randomA_sigmas.vector[seq.get,i],col=i,ylab="randomA_sigmas")
      }
      abline(h=0,col=(n.chains+1))
      
      dev.off()
      
    }
    
  }else{ #dumpConvergence=F
  
  out<-.C("GXE_mcmc", x=as.double(x), 
          n=as.integer(dim(x)[1]), as.double(y), 
          as.integer(N), as.integer(q), as.integer(num.haplo.id),
          as.integer(x.length), as.integer(h.length), 
          as.integer(freq.num), as.double(freq.new), 
          as.double(D), as.double(sigmas),as.double(beta), as.double(gamma), as.double(a), as.double(b),
          as.double(asig),as.double(bsig),  as.double(lambda), as.integer(num.it), 
          as.integer(burn.in), beta.out=as.double(beta.out), gamma.out=as.double(gamma.out),
          lambda.out=as.double(lambda.out),  sigmas.out=as.double(sigmas.out), randomA_sigmas.out=as.double(randomA_sigmas.out), freq.out=as.double(freq.out),D.out=as.double(D.out),as.integer(n.cov),up.xz=as.double(up.xz),as.integer(CC),BDIC.out=as.double(BDIC.out),
          as.double(randomA_sigmas),as.double(randomA_sigmas_a),as.double(randomA_sigmas_b),as.double(R_vec))
  
  
  
  
  BDIC.out<-out$BDIC.out
  BDIC_result=-2*mean(BDIC.out)+2*var(BDIC.out)
  
  
  
  #### what values we have in C
  cat("\n", "what values we have in C","\n")
  cat("Prior mean of beta.0 is ", out$beta0pm,"\n")
  cat("Prior variance of beta.0 is ", out$beta0pv,"\n")
  ####
  
  
  
  #### results from C program
  #up.xz<-matrix(out$up.xz,nrow=N,byrow=TRUE)
  #cat("head of up.xz is \n")
  #print(head(up.xz))
  
  #cat("The local address is",getwd(),"\n")
  #test<-list(x.t=x,n.t=dim(x)[1],y.t=y,N.t=N,num.haplo.id.t=num.haplo.id,
  #           x.length.t=x.length,h.length.t=h.length,freq.num.t=freq.num,
  #           freq.new.t=freq.new,D.t=D,sigmas.t=sigmas,beta.t=beta,a.t=a,b.t=b,
  #           asig.t=asig,bsig.t=bsig,lambda.t=lambda,num.it.t=num.it,burn.in.t=burn.in,
  #           beta.out.t=beta.out,lambda.out.t=lambda.out,sigmas.out.t=sigmas.out,
  #           freq.out.t=freq.out,D.out.t=D.out,n.cov.t=n.cov,up.xz.t=up.xz)
  #save(test,file=paste("testing_",SDD,"_file.RData",sep=""))
  #save.image()
  ####
  
  
  
  ###################################
  #####     start from here    ######
  ###################################
  
  ########OUTPUT########
  beta.out<-matrix(out$beta.out,nrow=num.it-burn.in, byrow=TRUE)
  gamma.out<-matrix(out$gamma.out,nrow=num.it-burn.in, byrow=TRUE)
  freq.out<-matrix(out$freq.out,nrow=num.it-burn.in, byrow=TRUE)
  sigmas.out<-out$sigmas.out
  randomA_sigmas.out<-out$randomA_sigmas.out
  lambda.out<-out$lambda.out
  D.out<-out$D.out
  
  
  #### output from LBLq_han.c
  cat("Some results from LBLq_han C function\n")
  cat("Mean sigma squared is ",mean(sigmas.out),"\n")
  cat("End Some results from LBLq_han C function\n")
  
  ####
  
  ci.beta<-numeric((x.length+1)*2)
  ci.gamma<-numeric(q*2)
  ci.lambda<-numeric(2)
  ci.D<-numeric(2)
  post.mean.beta<-numeric(x.length+1)
  post.mean.gamma<-numeric(q)
  ci.freq<-numeric((h.length)*2)
  post.mean.freq<-numeric(h.length)
  
  k<-1
  for (i in 1:(x.length+1))
  {
    ci.beta[k]<-quantile(beta.out[,i], probs=0.025)
    ci.beta[k+1]<-quantile(beta.out[,i], probs=0.975)
    k<-k+2
    post.mean.beta[i]<- mean(beta.out[,i])
  }
  
  k<-1
  for (i in 1:q)
  {
    ci.gamma[k]<-quantile(gamma.out[,i], probs=0.025) #vector
    ci.gamma[k+1]<-quantile(gamma.out[,i], probs=0.975)
    k<-k+2
    post.mean.gamma[i]<- mean(gamma.out[,i])
  }
  
  #OR<-exp(post.mean.beta)
  #OR<-round(OR,4)
  
  #### Take a look at what beta0 is
  cat("Mean of beta0 is ",post.mean.beta[1],"\n")
  cat("Lower 95 % CI of beta0 is ",ci.beta[1],"\n")
  cat("Upper 95 % CI of beta0 is ",ci.beta[2],"\n")
  ####
  
  
  ci.lambda<-c(quantile(out$lambda.out, probs=0.025),quantile(out$lambda.out, probs=0.975))
  ci.randomA_sigmas<-c(quantile(out$randomA_sigmas.out, probs=0.025),quantile(out$randomA_sigmas.out, probs=0.975))
  ci.sigmas<-c(quantile(out$sigmas.out, probs=0.025),quantile(out$sigmas.out, probs=0.975))
  ci.D<-c(quantile(out$D.out, probs=0.025),quantile(out$D.out, probs=0.975))
  
  
  ### check what CI of D is
  cat("Lower 95 % CI of D is ",ci.D[1],"\n")
  cat("Upper 95 % CI of D is ",ci.D[2],"\n")
  ######
  
  
  k<-1
  for (i in 1:(h.length))
  {
    ci.freq[k]<-quantile(freq.out[,i], probs=0.025)
    ci.freq[k+1]<-quantile(freq.out[,i], probs=0.975)
    k<-k+2
    post.mean.freq[i]<- mean(freq.out[,i])
  }
  
  #BF for beta
  prob.alt<-numeric(x.length+1)
  BF<-numeric(x.length+1)
  for (i in 1:(x.length+1))
  {
    prob.alt[i]<-length(beta.out[,i][abs(beta.out[,i]) > e])/(num.it-burn.in)
  }
  for (i in 1:(x.length+1))
  {
    prior.prob<-(b/(e+b))^a
    if (prob.alt[i]<1)
    {
      
      BF[i]<-(prob.alt[i]/(1-prob.alt[i]))/(prior.prob/(1-prior.prob))
      BF[i]<-round(as.numeric(BF[i]),2)
      if (BF[i]>100){
        BF[i]<-">100"
      }
    } else
      BF[i]<-">100"
  }
  
  
  #BF for gamma
  prob.alt<-numeric(q)
  BF.gamma<-numeric(q)
  for (i in 1:q)
  {
    prob.alt[i]<-length(gamma.out[,i][abs(gamma.out[,i]) > e])/(num.it-burn.in)
  }
  for (i in 1:q)
  {
      prior.prob<-0.9294654 # prior odds 13 in favor of Ha(not 0)
    if (prob.alt[i]<1)
    {
      
      BF.gamma[i]<-(prob.alt[i]/(1-prob.alt[i]))/(prior.prob/(1-prior.prob))
      BF.gamma[i]<-round(as.numeric(BF.gamma[i]),2)
      if (BF.gamma[i]>100){
          BF.gamma[i]<-">100"
      }
    } else
    BF.gamma[i]<-">100"
  }
  
  
  
  
  
  name1<-colnames(hdat[, -1,drop=F])
  #name2<-paste(name1[1],name1[-1],sep=".")
  pname<-c(baseline,name1)
  names(BF)<-pname
  names(post.mean.beta)<-pname
  #names(OR)<-pname
  if(q==1){
    names.gamma<-"PC1" #if R n by 1, name is PC1
  }else{
    names.gamma<-colnames(R)
  }
  names(BF.gamma)<-names.gamma
  names(post.mean.gamma)<-names.gamma

  
  ci.beta<-matrix(ci.beta,nrow=x.length+1, ncol=2, byrow=TRUE)
  #ci.OR<-exp(ci.beta)
  #ci.OR<-round(ci.OR,4)
  
  ci.beta<-data.frame(pname,ci.beta)
  colnames(ci.beta)<-c("Name", "Lower", "Upper")
  
  #ci.OR<-data.frame(pname,ci.OR)
  #colnames(ci.OR)<-c("Name", "Lower", "Upper")
  
  
  ci.gamma<-matrix(ci.gamma,nrow=q, ncol=2, byrow=TRUE)
  ci.gamma<-data.frame(names.gamma,ci.gamma)
  colnames(ci.gamma)<-c("Name", "Lower", "Upper")
  
  
  ci.freq<-matrix(ci.freq,nrow=(h.length), ncol=2, byrow=TRUE)
  ci.freq<-data.frame(c(colnames(hdat[, 2:h.length]),baseline), ci.freq)
  colnames(ci.freq)<-c("Hap", "Lower", "Upper")
  names(post.mean.freq)<-c(colnames(hdat[, 2:h.length]), baseline)
  
  
  BF<-c(BF,BF.gamma)
  PE<-c(post.mean.beta,post.mean.gamma)
  CI.PE<-rbind(ci.beta,ci.gamma)
  
  ans <- list(BF = BF, PE = PE, CI.PE = CI.PE,CI.lambda=ci.lambda, CI.sigmasA=ci.randomA_sigmas,CI.sigmas=ci.sigmas,freq=post.mean.freq)
  return(ans)
  rm(out)
  rm(BDIC.out)
  rm(freq.out) 
  rm(beta.out)
  rm(gamma.out)
  rm(D.out)
  rm(lambda.out)
  rm(sigmas.out)
  rm(randomA_sigmas.out)
 }
  
}
