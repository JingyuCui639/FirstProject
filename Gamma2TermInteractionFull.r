start.time <- Sys.time()
sink("GammaFullNewModel实验代码结果Interaction.txt")

#Finding the value of dispersion parameter p
#out<-tweedie.profile(Clinician$PracticeScore~Clinician$ebstatus+Clinician$monitoring+Clinician$implement_time, p.vec=seq(1.9, 5, .05), do.plot=T)
#Warnning message, but p.max=2

pterm<-2

#count the number of the objects, and give the number to "N"
N<-nrow(Clinician)

#initial the wards number "I" as 208 and doctors number "J" as 33
I<-208
J<-33

#construct indicator matrix B
 BI.mtx <- Matrix(0, N, I*J,sparse=TRUE)
 BJ.mtx <- Matrix(0, N, I*J,sparse=TRUE)
 
 for (k in 1:N)
 {
   BI.mtx[k,(Clinician$Subject[k]-1)*J+Clinician$Supervison[k]]<-1
   BJ.mtx[k,(Clinician$Supervison[k]-1)*I+Clinician$Subject[k]]<-1
 }  
 
 #Creating Covariates Matrix
 X<-cbind(1,Clinician$ebstatus,Clinician$monitoring,Clinician$time1,Clinician$time2,Clinician$time3,Clinician$ebstatus*Clinician$monitoring,
 Clinician$ebstatus*Clinician$time1,Clinician$ebstatus*Clinician$time2,Clinician$ebstatus*Clinician$time3,Clinician$monitoring*Clinician$time1,
 Clinician$monitoring*Clinician$time2,Clinician$monitoring*Clinician$time3)
 colnames(X)<-c("constant","ebstatus","monitoring","time1","time2","time3","ebstatus*monitoring",
 "ebstatus*time1","ebstatus*time2","ebstatus*time3","monitoring*time1","monitoring*time2","monitoring*time3")
 
 #The number of regression parameters P
 #P<-ncol(X)
 
#Creating y--response vextor
y<-5-Clinician$PracticeScore
 
#Calculate the best estimates of beta using GLM
fit <- glm(y ~ X[,2]+X[,3]+X[,4]+X[,5]+X[,6]+X[,7]+X[,8]+X[,9]+X[,10]+X[,11]+X[,12]+X[,13] , family = Gamma)
summary(fit)
 print("summary(fit)$coefficients")
 summary(fit)$coefficients
 
#-----Set initial values---------------------------------

print("Initial values")
betaest<-matrix(c(summary(fit)$coefficients[,1]))
print("betaest")
print(betaest)

# Calculate the unconditioal mean for y
  mu.uncondi<-exp(X%*%betaest) 
#  print("mu.uncondi")
#  print(mu.uncondi)
  diag_mu<-Matrix((diag(as.vector(mu.uncondi))),sparse=TRUE)
  diag_mu.1_pterm<-Matrix(diag(as.vector(mu.uncondi^(1-pterm))),sparse=TRUE)

#Overall response average
ybar<-mean(y)

#Initial value of U---I*1 matrix
U<-matrix(0, nrow=I)

 for (i in 1:I){
     U[i]<-mean(y[rowSums(BI.mtx[,c(((i-1)*J+1):(i*J))])=="1"])/ybar
	 }

#Initial value of V----J*1 matrix
V<-matrix(0, nrow=J)

for (i in 1:J){
    V[i]<-mean(y[rowSums(BJ.mtx[,c(((i-1)*I+1):(i*I))])=="1"])/ybar
	}

#Initial value of U*hat---------(I*J)*1 matrix
  
ustar.hat<-matrix(0,nrow=(I*J),ncol=1)	
for (i in 1:(I*J)){
    if(sum(BI.mtx[,i]!=0)){
	    ustar.hat[i]<-mean(y[BI.mtx[,i]=="1"])/ybar
	    }
	}
	
#Calculating conditional mean
ustar.data<-BI.mtx%*%ustar.hat
mu.condi<-mu.uncondi*ustar.data
#diag_mu.condi<-Matrix(diag(as.vector(mu.condi)),sparse=TRUE)
#diag_mu.condi.1_pterm<-Matrix(diag(as.vector(mu.condi^(1-pterm))),sparse=TRUE)
	
#Set the initial value of sigmasq 

sigmasq<-mean((U-1)^2)
print("sigmasq")
print(sigmasq)

#set the initial value of tausq
tausq<-mean((V-1)^2)
print("tausq")
print(tausq)

#Set the initial value of rousq

 rousq<-mean(((y-mu.condi)^2)/(mu.uncondi^pterm))

 print("rousq") 
 print(rousq)
 
 betaest.list<-list()
 sigmasq.list<-list()
 tausq.list<-list()
 rousq.list<-list()
 se.beta.list<-list()
 diff.list<-list()
 
################################################################
#Start iterations
Times<-0
diff<-1
while(diff>=10^(-5)){
 
 Times<-Times+1
  print("iterition time:")
  print(Times)
  
  betaest_old<-betaest
  sigmasq_old<-sigmasq
  tausq_old<-tausq
  rousq_old<-rousq

# --------------------Estimating new U* hat------------------------------------ 

#Calculating Var(U*)
#constructing layer1 of V(U*)
	layer1.block<-as.matrix(tausq*Diagonal(J))
	list1<-mget(rep("layer1.block",I))
	layer1.row<-do.call("cbind",list1)
	list2<-mget(rep("layer1.row",I))
	layer1<-Matrix(do.call("rbind",list2))
#constructing layer2 of V(U*)
    layer2.block<-sigmasq*matrix(1,nrow=J,ncol=J)
	layer2<-bdiag(replicate(I,layer2.block,simplify=F))
#constructing layer3 of V(U*)
    layer3<-sigmasq*tausq*Diagonal(I*J)
#Final result of V(U*)
    var.ustar<-layer1+layer2+layer3

#write.csv(var.ustar, file="I:/STAT_UNBF/Thesis/R Program/Pupil score/var.ustar.csv")

#create covariance of U* and Y---cov.ustary
cov.ustary<-var.ustar%*%t(BI.mtx)%*%diag_mu

# create matrix var(Y)
var.Y<-rousq*(diag_mu^pterm)+diag_mu%*%BI.mtx%*%cov.ustary
var.Y.inv<-solve(var.Y)

#Calculating U* hat
ustar.hat<-as.matrix(matrix(1,nrow=(I*J))+cov.ustary%*%var.Y.inv%*%(matrix(y)-mu.uncondi))
#print("U* hat")
#print(ustar.hat)

#Calculate Var(U* hat)
var.ustar.hat<-cov.ustary%*%var.Y.inv%*%t(cov.ustary)

#calculate conditional mean 
ustar.data<-BI.mtx%*%ustar.hat
mu.condi<-mu.uncondi*ustar.data
#diag_mu.condi<-Matrix(diag(as.vector(mu.condi)),sparse=TRUE)

#-----------------------------------Estimate new beta----------------------------
#Creating Psi beta function---calculating in matrix multiplication
Psi.beta<-as.matrix(t(X)%*%diag_mu.1_pterm%*%(matrix(y)-mu.condi)/rousq)
print("Psi.beta")
print(Psi.beta)

#Creating S(beta)
#S.beta<-as.matrix(-t(X)%*%diag_mu.1_pterm%*%diag_mu.condi%*%X/rousq)
V.beta<-as.matrix(t(X)%*%diag_mu%*%var.Y.inv%*%diag_mu%*%X)
print("V(beta)")
print(V.beta)
S.beta<--V.beta
print("S.beta")
print(S.beta)

S.betainv<-solve(S.beta)
print("S.betainv")
print(S.betainv)

#Calculate the covariance matrix of beta
var.beta<-solve(V.beta)
print("var.beta")
print(var.beta)
# Calculating standard error of betaest
std.error<-sqrt(diag(var.beta))
se.beta.list[[Times]]<-std.error
print("s.e. of old betaest")
print(std.error)

zvalue<- matrix(betaest)/matrix(std.error)
print("z- value of old betaest (different from 0):")
print(zvalue)

#Finding the new beta estimate
betaest<-betaest_old-S.betainv%*%Psi.beta
betaest.list[[Times]]<-betaest

print("new betaest:")
print(betaest)

#difference_beta<-betaest-betaest_old
#print("difference between new and old beta:")
#print(difference_beta)

# Calculate the new unconditioal mean for y
   mu.uncondi<-exp(X%*%betaest) 
#  print("mu.uncondi")
#  print(mu.uncondi)
  diag_mu<-Matrix((diag(as.vector(mu.uncondi))),sparse=TRUE)
  diag_mu.1_pterm<-Matrix(diag(as.vector(mu.uncondi^(1-pterm))),sparse=TRUE)
  
#calculate conditional mean 
#ustar.data<-BI.mtx%*%ustar.hat
mu.condi<-mu.uncondi*ustar.data
#diag_mu.condi<-Matrix(diag(as.vector(mu.condi)),sparse=TRUE)

#-------------------------Estimat U  hat new-------------------------
#Calculating cov(U,U*)

vector.sigmasq<-sigmasq*matrix(1,ncol=J)
cov.uustar<-bdiag(replicate(I,vector.sigmasq,simplify=F))

#write.csv(cov.uustar, file="I:/STAT_UNBF/Thesis/R Program/Pupil score/cov.uustar.csv")

cov.uy<-cov.uustar%*%t(BI.mtx)%*%diag_mu
uhat<-as.matrix(matrix(1,nrow=I)+cov.uy%*%var.Y.inv%*%(matrix(y)-mu.uncondi))
#print("U hat")
#print(uhat)

#---------------------------Estimate V hat new--------------------------

block.tausq<-as.matrix(tausq*Diagonal(J))
list<-mget(rep("block.tausq",I))
cov.vustar<-do.call("cbind",list)

#write.csv(cov.vustar, file="I:/STAT_UNBF/Thesis/R Program/cov.vustar.csv")
cov.vy<-cov.vustar%*%t(BI.mtx)%*%diag_mu
vhat<-as.matrix(matrix(1,nrow=J)+cov.vy%*%var.Y.inv%*%(matrix(y)-mu.uncondi))
#print("V hat")
#print(vhat)

#-------------------------Estimate new sigmasq--------------------------------------------

sumI<-0
for (i in 1:I){
    nrow<-c(((i-1)*J+1):(i*J))
    #BlockI<-var.ustar.hat[nrow,nrow]
    sumI<-sumI+sum(var.ustar.hat[nrow,nrow])
}
sumI<-sumI-sum(diag(var.ustar.hat))
Ci<-sigmasq-(sumI/(I*J*(J-1)))

sigmasq<-mean((uhat-1)^2)+Ci
sigmasq.list[[Times]]<-sigmasq
#difference_sigmasq<-sigmasq-sigmasq_old
print("new sigmasq")
print(sigmasq)
#print("Difference between new and old sigma:")
#print(difference_sigmasq)

#----------------Estimating new tausq----------------------------------

sumJ<-0
for (i in 1:(I-1)){
   nrow<-c(1:(i*J))
   ncol<-c(((I-i)*J+1):(I*J))
   sumJ<-sumJ+sum(diag(var.ustar.hat[nrow,ncol]))
}
meanJ<-sumJ/((I*(I-1)*J)/2)
Cj<-tausq-meanJ

tausq<-mean((vhat-1)^2)+Cj
tausq.list[[Times]]<-tausq
#difference_tausq<-tausq-tausq_old
print("new tausq")
print(tausq)
#print("difference between new and old tausq:")
#print(difference_tausq)

#-----------------------Estimating new rousq --------------------------------

part1<-((y-mu.condi)^2)/(mu.uncondi^pterm)
D<-Diagonal(x=diag(var.ustar.hat))
temp<-BI.mtx%*%D%*%t(BI.mtx)
diag_Cstar.data<-(sigmasq+tausq+sigmasq*tausq)-diag(temp)
part2<-mu.uncondi^(2-pterm)*matrix(diag_Cstar.data)
rousq<-mean(part1+part2)

 rousq.list[[Times]]<-rousq
 
 print("new rousq:")
 print(rousq)
 #print("difference between old and new rousq:")
 #print(rousq_old-rousq)



diff<-abs(sigmasq-sigmasq_old)+abs(tausq-tausq_old)+abs(rousq-rousq_old)+sum(abs(betaest-betaest_old))
diff.list[[Times]]<-diff
print("difference")
print(diff)

}
print("#The end of iteration")
##########################################################################

#Calculating Var(U*)
#constructing layer1 of V(U*)
	layer1.block<-as.matrix(tausq*Diagonal(J))
	list1<-mget(rep("layer1.block",I))
	layer1.row<-do.call("cbind",list1)
	list2<-mget(rep("layer1.row",I))
	layer1<-Matrix(do.call("rbind",list2))
#constructing layer2 of V(U*)
    layer2.block<-sigmasq*matrix(1,nrow=J,ncol=J)
	layer2<-bdiag(replicate(I,layer2.block,simplify=F))
#constructing layer3 of V(U*)
    layer3<-sigmasq*tausq*Diagonal(I*J)
#Final result of V(U*)
    var.ustar<-layer1+layer2+layer3
		
#create covariance of U* and Y---cov.ustary
cov.ustary<-var.ustar%*%t(BI.mtx)%*%diag_mu

# create matrix var(Y)
var.Y<-rousq*(diag_mu^pterm)+diag_mu%*%BI.mtx%*%cov.ustary
var.Y.inv<-solve(var.Y)
#--------------------------------------------------------------------
#Creating S(beta)
#S.beta<-as.matrix(-t(X)%*%diag_mu.1_pterm%*%diag_mu.condi%*%X/rousq)
V.beta<-as.matrix(t(X)%*%diag_mu%*%var.Y.inv%*%diag_mu%*%X)
print("V.beta")
print(V.beta)

#Calculate the covariance matrix of beta
var.beta<-solve(V.beta)
print("var.beta")
print(var.beta)
# Calculating standard error of  new betaest
std.error<-sqrt(diag(var.beta))
se.beta.list[[Times+1]]<-std.error
print("s.e. of new estimated beta")
print(std.error)

zvalue<- betaest/matrix(std.error)
print("z- value of new estimated beta (different from 0):")
print(zvalue)

pvalues<-2*pnorm(-abs(zvalue))
print("pvalues for betaest")
print(pvalues)

print("Total number of iteration:")	
print(Times)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
#betaest.list
#sigmasq.list
#tausq.list
sink()


Ite.Num<-Times 
#Drawing the picture of how the estimates of constant changes
beta0<-matrix(0,nrow = Ite.Num)
for (i in 1:Ite.Num)
{
   beta0[i,]<-betaest.list[[i]][1,]
 }
 plot(beta0,main="The changing of beta0",ylab="Estimates of beta0",type="l",col="blue")
 
#Drawing the picture of how the estimates of beta1 changes
beta1<-matrix(0,nrow = Ite.Num)
for (i in 1:Ite.Num)
{
   beta1[i,]<-betaest.list[[i]][2,]
 }
 plot(beta1,main="The changing of beta1",ylab="Estimates",type="l",col="red")
 
 #Drawing the picture of how the estimates of beta2 changes
beta2<-matrix(0,nrow = Ite.Num)
for (i in 1:Ite.Num)
{
   beta2[i,]<-betaest.list[[i]][3,]
 }
 plot(beta2,main="The changing of beta2",ylab="Estimates",type="l",col="green")

 #Drawing the picture of how the estimates of beta3 changes
beta3<-matrix(0,nrow = Ite.Num)
for (i in 1:Ite.Num)
{
   beta3[i,]<-betaest.list[[i]][4,]
 }
 plot(beta3,main="The changing of beta3",ylab="Estimates",type="l",col="black")

 plot(as.matrix(sigmasq.list),main="The changing of sigmasq",ylab="Estimates",type="l",col="brown")
 plot(as.matrix(tausq.list),main="The changing of tausq",ylab="Estimates",type="l",col="orange")
 plot(as.matrix(rousq.list),main="The changing of rousq",ylab="Estimates",type="l",col="black")

 #Draw the picture of s.e. of beta0
 se.beta0<-rep(0,Ite.Num+1)
 for (i in 1:Ite.Num+1){
     se.beta0[i]<-se.beta.list[[i]][1]
}
plot(se.beta0, main="changing of s.e. of beta0", ylab="estimates", type="l",col="blue")

#Draw the picture of s.e. of beta1
 se.beta1<-rep(0,Ite.Num+1)
 for (i in 1:Ite.Num+1){
     se.beta1[i]<-se.beta.list[[i]][2]
}
plot(se.beta1, main="changing of s.e. of beta1", ylab="estimates", type="l",col="red")

#Draw the picture of s.e. of beta2
 se.beta2<-rep(0,Ite.Num+1)
 for (i in 1:Ite.Num+1){
     se.beta2[i]<-se.beta.list[[i]][3]
}
plot(se.beta2, main="changing of s.e. of beta2", ylab="estimates", type="l",col="brown")