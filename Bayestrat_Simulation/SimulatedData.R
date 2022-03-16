# Example of simulation for population structure setting 1. 
#This is the power setting where the first 25 SNPs are associated with the phenotype.

#### simulate ancestry and population allele frequencies ####
set.seed(12)
n.all<-10000 #number of null and candidate SNPs
F_st<-0.01

p_anc<-runif(n.all,0.05,0.95)

p_sub<-matrix(0,nrow=n.all,ncol=4)
for(i in 1:n.all){
	p_sub[i,]<-rbeta(4,(1-F_st)*p_anc[i]/F_st,(1-F_st)*(1-p_anc[i])/F_st)
}

for(j in 1:4){
	for(i in 1:n.all){
		if(p_sub[i,j]>0.5){
			p_sub[i,j]<-1-p_sub[i,j] #p mimor allele frequency
		}
	}
}


#### set effect sizes ####
SNP.ind<-seq(1,100) #reserve the first 100 SNPs for testing
n.SNP<-length(SNP.ind)
n<-250 #250 individuals in each of the 4 populations, total sample size=250*4=1000

mu<-c(0,9,6,3) #mean shift parameter in simulation structure 1, notice the order of the 4 populations were changed.
SNP.eff.id<-c(1:25) #the first 25 PCs are associated
c<-0.15 #this is power, for type I error set c=0
beta<-c*abs(log(apply(p_sub[SNP.ind[SNP.eff.id],],1,mean)))
beta #SNP effect size

#### simulate phenotype ####
pop1<-matrix(0,nrow=n,ncol=1+n.SNP)
pop2<-matrix(0,nrow=n,ncol=1+n.SNP)
pop3<-matrix(0,nrow=n,ncol=1+n.SNP)
pop4<-matrix(0,nrow=n,ncol=1+n.SNP)

for(i in 1:n){
	X<-sapply(SNP.ind,function(x){sum(sample(c(0,1),2,prob=c(1-p_sub[x,1],p_sub[x,1]),replace=T))},simplify=T)
	pop1[i,1]<-mu[1]+X[SNP.eff.id]%*%beta+rnorm(1)  #Y=mu+snp*beta+epsilon, right now SNP.eff.id has to be within SNP.id
	pop1[i,-1]<-X
	
	X<-sapply(SNP.ind,function(x){sum(sample(c(0,1),2,prob=c(1-p_sub[x,2],p_sub[x,2]),replace=T))},simplify=T)
	pop2[i,1]<-mu[2]+X[SNP.eff.id]%*%beta+rnorm(1)  #Y=mu+snp*beta+epsilon
	pop2[i,-1]<-X
	
	X<-sapply(SNP.ind,function(x){sum(sample(c(0,1),2,prob=c(1-p_sub[x,3],p_sub[x,3]),replace=T))},simplify=T)
	pop3[i,1]<-mu[3]+X[SNP.eff.id]%*%beta+rnorm(1)  #Y=mu+snp*beta+epsilon
	pop3[i,-1]<-X
	
	X<-sapply(SNP.ind,function(x){sum(sample(c(0,1),2,prob=c(1-p_sub[x,4],p_sub[x,4]),replace=T))},simplify=T)
	pop4[i,1]<-mu[4]+X[SNP.eff.id]%*%beta+rnorm(1)  #Y=mu+snp*beta+epsilon
	pop4[i,-1]<-X
	
}

dat<-rbind(pop1,pop2,pop3,pop4)
dat<-as.data.frame(dat)
colnames(dat)[1]<-"Y"
save(dat,file="SimulatedDataStruc1Power.rda")
	




