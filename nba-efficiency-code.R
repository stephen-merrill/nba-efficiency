# load and preprocess data, available at https://www.kaggle.com/dansbecker/nba-shot-logs

shots <- read.csv('C:/Users/Stephen/Google Drive/651/shot_logs.csv')
shots <- shots[shots$TOUCH_TIME > 0,]
shots[is.na(shots)] <- 0 

# this is a shorter analysis with only three players for computational benefits
players <- c("gordon hayward","kobe bryant","jimmer dredette")
shots <- shots[shots$player_name %in% players,]
n.players <- length(unique(shots$player_name))
player.names <- unique(shots$player_name)
library(dplyr)
n.player.shots <- count(shots,player_name)$n
#sort alphabetically
shots <- shots[order(shots$player_name),]
player.names <- sort(player.names)
#change H to 1 and A to 0
shots$LOCATION <- as.numeric(shots$LOCATION)
shots$LOCATION <- shots$LOCATION-1

####################################################################
############                 Model              ####################
####################################################################
# Bayesian Hierarchical Latent Variable Probit Model

library(msm)
library(MCMCpack)

#this is a shorter mcmc run, still has nice convergence
length<-400
burn<-100
y<-shots$PTS

# variables of interest:
# LOCATION
# SHOT_CLOCK
# DRIBBLES
# TOUCH_TIME
# CLOSE_DEF_DIST

X.list <- list()
for(i in 1:n.players) {
  X.list[[i]] <- shots[shots$player_name==unique(shots$player_name)[i],c("LOCATION","SHOT_CLOCK","DRIBBLES","TOUCH_TIME","CLOSE_DEF_DIST")]
}

#X.list has n.players # of matricies, each with a different n # of rows and 5 columns
#n x 5 x 146

#setting up for mcmc params
beta<-matrix(0,ncol(X.list[[1]]),nrow=n.players)
beta.save<-array(dim = c(n.players,ncol(X.list[[1]]),(burn+length)))
beta.save[,,1]<-beta

mu<-matrix(0,ncol=ncol(X.list[[1]]),nrow=(burn+length))
mu[1,]<-rep(0,ncol(X.list[[1]]))

m<-rep(0,ncol(X.list[[1]]))
V<-1*diag(ncol(X.list[[1]]))

w<-1+n.players
I<-diag(ncol(X.list[[1]]))
Sigma<-array(dim = c(ncol(X.list[[1]]),ncol(X.list[[1]]),(burn+length)))
Sigma[,,1]<-I

Z.list <- list()
for(i in 1:n.players) {
  Z.list[[i]] <- matrix(0,nrow=(length+burn),ncol=n.player.shots[i])
}
#length+burn x n.player.shots x n.players

#gamma cut points are fixed logically
gamma<-c(2,3)

devs <- matrix(0,ncol=n.players,nrow=length+burn)

#run mcmc
for(d in 2:(length+burn)){
  # update beta
  for (k in 1:n.players){
    sigstar <- solve( solve(Sigma[,,(d-1)]) + t(X.list[[k]])%*%as.matrix(X.list[[k]]) )
    mustar <- sigstar %*% (solve(Sigma[,,(d-1)])%*%mu[(d-1),] + t(X.list[[k]])%*%Z.list[[k]][(d-1),])
    beta.save[k,,d] <- mvrnorm(1,mustar,sigstar)
  }
  
  # update mu
  Vstar <- solve( solve(V) +  solve(Sigma[,,(d-1)]) )
  mstar <- Vstar %*% (solve(V)%*%m + solve(Sigma[,,(d-1)])%*%colMeans(beta.save[,,d]))
  mu[d,] <- mvrnorm(1,mstar,Vstar)
  
  # update Sigma
  wstar<-w+n.players
  S.mu<-matrix(0,nrow = ncol(X.list[[1]]),ncol = ncol(X.list[[1]]))
  for(k in 1:n.players){
    S.mu<-S.mu+(beta.save[k,,d]-mu[d,])%*%t((beta.save[k,,d]-mu[d,]))
  }
  Istar<-solve(I+S.mu)
  Sigma[,,d] <- riwish(wstar,Istar)
  
  # update Z
  index <- 1
  for(k in 1:n.players){ 
    for (i in 1:n.player.shots[k]){
      Z.list[[k]][d,i] <-rtnorm(1,t(X.list[[k]][i,])%*%beta.save[k,,d],1,gamma[2],Inf)
      if(y[index]==2)
      {
        Z.list[[k]][d,i] <-rtnorm(1,t(X.list[[k]][i,])%*%beta.save[k,,d],1,gamma[1],gamma[2])
      }
      if(y[index]==0)
      {
        Z.list[[k]][d,i] <-rtnorm(1,t(X.list[[k]][i,])%*%beta.save[k,,d],1,-Inf,gamma[1])
      }
      index <- index+1
    }
  }
  index <- 1
  for(k in 1:n.players) {
    p <- hist(Z.list[[k]][d,],plot = F,breaks=c(-Inf,2,3,Inf))$count/n.player.shots[[k]]
    devs[d,k] <- -2*loglike(y[index:(index+n.player.shots[k]-1)],p)
    index <- index + n.player.shots[k]
  }
}

loglike <- function(y,p) {
  counts <- c(sum(y[which(y==0)]),sum(y[which(y==2)]),sum(y[which(y==3)]))
  sum(dmultinom(counts,prob=p,log=TRUE))
}

#Model fit:

#attempt at dic
dbar <- mean(apply(devs[-c(1:burn),],1,sum))
dthetabar <- numeric()
index <- 1
for(k in 1:n.players) {
  p <- hist(apply(Z.list[[k]],2,mean),plot = F,breaks=c(-Inf,2,3,Inf))$count/n.player.shots[[k]]
  dthetabar[k] <- -2*loglike(y[index:(index+n.player.shots[k]-1)],p)
  index <- index + n.player.shots[k]
}
psubd <- dbar-sum(dthetabar)

#Chi Square goodness of fit
bincnts<-0
nbin<-floor(n.player.shots^0.4)
quants<-matrix(0,nrow = n.players,ncol = 4)
index <-1
for(k in 1:n.players){
  quants[k,]<-c(0,as.numeric(cumsum(table(y[index:(index+n.player.shots[k]-1)])/n.player.shots[k])))
  index <- index + n.player.shots[k]
}
BX2<-matrix(0,nrow = (length+burn),ncol = n.players)

for(d in 1:(length+burn)){
  index <- 1
  for(k in 1:n.players){
    cdf<-numeric(n.player.shots[k])
    for(i in 1:n.player.shots[k]){
      cdf[i]<-pnorm(0,t(X.list[[k]][i,])%*%beta.save[k,,(d+burn)],1)
      if(y[index]==2){
        cdf[i]<-pnorm(1,t(X.list[[k]][i,])%*%beta.save[k,,(d+burn)],1)
      }
      if(y[index]==3){
        cdf[i]<-1
      }
      index <- index + 1
    }
    bincnts<-hist(cdf,plot=F,breaks=quants[k,],na.rm=T)$count
    n<-sum(bincnts)
    BX2[d,k]<-sum((bincnts-n/nbin[k])^2/(n/nbin[k]))
  }
}

#Results:

#post pred
z.pred<-matrix(0,nrow = length,ncol = n.players)

X.samp<-array(dim = c(n.players,5,length))
for(j in 1:5){
  for(k in 1:n.players){
    X.samp[k,j,]<-sample(X.list[[k]][,j],length,replace = TRUE)
  }
}

for(k in 1:n.players){
  index<-sample((burn+1):M,length,replace = T)
  for(d in 1:length){
    z.pred[d,k]<-rnorm(1,t(X.samp[k,,d])%*%beta.save[k,,index[d]],1)
  }
}

prob.pred<-matrix(0,nrow = n.players,ncol = 3)

for(k in 1:n.players){
  prob.pred[k,]<-hist(z.pred[,k],plot = F,breaks=c(-Inf,2,3,Inf))$count/length
}

pred.avg.score<-numeric(9)
values<-matrix(c(0,2,3),ncol = 1)
pred.avg.score<-(prob.pred%*%values)

player.names <- as.character(player.names)
player.names[4] <- "jimmer fredette"

#plot betas trace
par(mfrow=c(4,1))
for(i in 1:n.players) {
  plot(beta.save[i,1,-c(1:burn)],type='l',main=player.names[i],ylab="LOC")
  plot(beta.save[i,2,-c(1:burn)],type='l',ylab="CLOCK")
  plot(beta.save[i,3,-c(1:burn)],type='l',ylab="DRIBS")
  plot(beta.save[i,4,-c(1:burn)],type='l',ylab="TIME")
  plot(beta.save[i,5,-c(1:burn)],type='l',ylab="DIST")
  #plot(beta.save[i,6,],type='l')
}

#plot beta densities
par(mfrow=c(1,1))
plot(density(beta.save[1,1,-c(1:burn)]),main="Location, Effect of Playing at Home",lwd=2,ylim=c(0,3.5),col="brown")
for(i in 2:n.players) {
  lines(density(beta.save[i,1,-c(1:burn)]),col=i,lwd=2)
}
legend("topleft",legend=player.names,lty=1,col=c("brown",2:n.players))

plot(density(beta.save[1,2,-c(1:burn)]),main="Shot Clock Time Remaining",lwd=2,col="brown",xlim=c(-0.03,0.1),ylim=c(0,80))
for(i in 2:n.players) {
  lines(density(beta.save[i,2,-c(1:burn)]),col=i,lwd=2)
}
legend("topleft",legend=player.names,lty=1,col=c("brown",2:n.players))

plot(density(beta.save[1,3,-c(1:burn)]),main="Number of Dribbles",lwd=2,ylim=c(0,20),col="brown")
for(i in 2:n.players) {
  lines(density(beta.save[i,3,-c(1:burn)]),col=i,lwd=2)
}
legend("topleft",legend=player.names,lty=1,col=c("brown",2:n.players))

plot(density(beta.save[1,4,-c(1:burn)]),main="Touch Time",ylim=c(0,20),lwd=2,col="brown")
for(i in 2:n.players) {
  lines(density(beta.save[i,4,-c(1:burn)]),col=i,lwd=2)
}
legend("topleft",legend=player.names,lty=1,col=c("brown",2:n.players))

plot(density(beta.save[1,5,-c(1:burn)]),xlim=c(-0.04,0.2),main="Closest Defender Distance",ylim=c(0,35),lwd=2,col="brown")
for(i in 2:n.players) {
  lines(density(beta.save[i,5,-c(1:burn)]),col=i,lwd=2)
}
legend("topleft",legend=player.names,cex=0.8,lty=1,col=c("brown",2:n.players))

