########################################
## Code for Dutch boys example
########################################

rm(list = ls())
detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}
detach_package("mice")
library(mice)
library(rjags)
library(mvtnorm)
summary(boys$age)
dim(boys)
boys<-boys[boys[,"age"]>=1 & boys[,"age"]<=18,]

boys[1:20,]
dim(boys) 

summary(boys$age)
boys$logbmi <- log(boys$bmi)
boys$loghgt <- log(boys$hgt)
boys$logwgt <- log(boys$wgt)
boys$city <- as.numeric(boys$reg=="city")
boys<-(boys[,c("hgt", "wgt", "logbmi", "loghgt", "logwgt", "city", "age")])

boxplot(boys$logbmi~boys$reg)

library(ggvenn)
ggvenn(as.data.frame(is.na(boys[,c("logbmi","loghgt","logwgt","city")])), 
fill_color=c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"))

# Frequentist analysis:
plot(logbmi~age,data=boys[boys[,"age"]>1,], 
    col=boys[boys[,"age"]>1,"city"]+1, pch=20)
mod1 <- lm(logbmi~city*age + I(age^2) ,
    data=boys[boys[,"age"]>1,])
freq_complete <- round(c(coef(mod1)["city"], confint(mod1)["city",]),3)
freq_complete
summary(mod1)


##############################
nMCMC <- 2000
thedata <- boys[,c("logbmi", "hgt", "wgt", "city", "age")]
dim(thedata)




# Univariate analysis:
# model 1: a model with y is input
jags_univ <- "model {
# Priors
beta0 ~ dnorm(0,1)
beta1 ~ dnorm(0,1)
beta2 ~ dnorm(0,1)
beta3 ~ dnorm(0,1)
beta4 ~ dnorm(0,1)
sigma ~ dexp(1)

# Model	
for(i in 1:N){
	y[i]     ~ dnorm(beta0 + beta1*x1[i] + 
                        beta2*x2[i]+ beta3*x1[i]*x2[i] + beta4*(x2[i]^2),
	 1/(sigma_squared))}

# Output
sigma_squared <- sigma^2
}"




# Bivariate analysis:
# model 2: a model with z1 and z2 as input
jags_biv<- "model {
# Priors
alpha0 ~ dnorm(0,1)
alpha1 ~ dnorm(0,1)
alpha2 ~ dnorm(0,1)
alpha3 ~ dnorm(0,1)
alpha4 ~ dnorm(0,1)

gamma0 ~ dnorm(0,1)
gamma1 ~ dnorm(0,1)
gamma2 ~ dnorm(0,1)
gamma3 ~ dnorm(0,1)
gamma4 ~ dnorm(0,1)

# Constructing the covariance matrix and the corresponding precision matrix.
    prec[1:2,1:2] <- inverse(cov[,])
    cov[1,1] <- sigma[1] * sigma[1]
    cov[1,2] <- sigma[1] * sigma[2] * rho
    cov[2,1] <- sigma[1] * sigma[2] * rho
    cov[2,2] <- sigma[2] * sigma[2]
    
# Flat priors on all parameters.
    sigma[1] ~ dexp(1) 
    sigma[2] ~ dexp(1) 
    rho ~ dunif(-1, 1)
    pi  ~ dunif(0, 1)
# Model	
for(i in 1:N){
  x1[i]  ~ dbin(pi,1);
  mu[i,1] = alpha0 + alpha1*x1[i] + 
  				alpha2*x2[i] + alpha3*x1[i]*x2[i] + alpha4*(x2[i]^2) ;
  mu[i,2] = gamma0 + gamma1*x1[i] + 
  				gamma2*x2[i] + gamma3*x1[i]*x2[i] + gamma4*(x2[i]^2);
	#z[i,1:2] ~ dmnorm(mu[i,1:2], prec[1:2,1:2]);
	
	
	z[i,1] ~ dmnorm(mu[i,1], pow(sigma[1],-2) );
	z[i,2] ~ dmnorm(mu[i,2] + (sigma[2]/sigma[1])*rho*(z[i,1]-mu[i,1]), 
	                1/((1-(rho^2))*(sigma[2]^2)) );
	
}

}"


############################################################
############################################################
### ### With univariate ### ###
start_time <- Sys.time()
# model with y as input
thedata_cc<-na.omit(thedata)
jags.m <- jags.model(textConnection(jags_univ), 
                     data = list(y = thedata_cc[,"logbmi"], 
                                 x1 = thedata_cc[,"city"],
                                 x2 = thedata_cc[,"age"],
                                 N = dim(thedata_cc)[1]))

# this is our estimate
mAsamples <-(coda.samples(jags.m, c("beta1", "beta3"),
    n.iter = nMCMC, n.burnin=1000))
theta_samples <- apply(mAsamples[[1]],1, function(q) {
q["beta1"] + q["beta3"]*mean(thedata[,"age"])})  

univariate_complete_est <- round(quantile(
        theta_samples, c(0.5,0.025,0.975)),3)
end_time <- Sys.time()
mAtime <- end_time-start_time
univariate_complete_time <-mAtime
univariate_complete_est
univariate_complete_time
############################################################
############################################################

####################
create_theta_gcomp <- function(q){

#	N <- dim(thedata)[1]
	N <- 1
	x1 = rep(0,S*N)
	x2 = sample(thedata[,"age"], S*N, replace=TRUE)
	
		mu1 = q["alpha0"] + q["alpha1"]*x1 + 
                q["alpha2"]*x2 + q["alpha3"]*x1*x2 + q["alpha4"]*(x2^2);
		mu2 = q["gamma0"] + q["gamma1"]*x1 + 
                q["gamma2"]*x2 + q["gamma3"]*x1*x2 + q["gamma4"]*(x2^2);

		Sigma = matrix(c(q["sigma[1]"]^2, 
		q["sigma[1]"]*q["sigma[2]"]*q["rho"],
		q["sigma[1]"]*q["sigma[2]"]*q["rho"],
		q["sigma[2]"]^2),2,2) 
	
	
logBMIstar0 <- apply(cbind(1:(S*N)),1, function(i){
			zstar_i <- rmvnorm(1,cbind(mu1[i], mu2[i]), sigma= Sigma)
			logBMIstar <- (zstar_i[2]) - log((exp(zstar_i[1])/100)^2)
			return(logBMIstar)})

# Y* ~ Y|X1=0,X2=samplex2
			
	x1 = rep(1,S*N)
	
		mu1 = q["alpha0"] + q["alpha1"]*x1 + 
                q["alpha2"]*x2 + q["alpha3"]*x1*x2 + q["alpha4"]*(x2^2);
		mu2 = q["gamma0"] + q["gamma1"]*x1 + 
                q["gamma2"]*x2 + q["gamma3"]*x1*x2 + q["gamma4"]*(x2^2);

	
logBMIstar1 <- apply(cbind(1:(S*N)),1, function(i){
			zstar_i <- rmvnorm(1,cbind(mu1[i], mu2[i]), sigma= Sigma)
			logBMIstar <- (zstar_i[2]) - log((exp(zstar_i[1])/100)^2)
			return(logBMIstar)})	
# Y* ~ Y|X1=1,X2=samplex2			
			
	thetasample <- mean(logBMIstar1)  -	mean(logBMIstar0)	
	
	return(thetasample)}
########################################


############################################################
### ### With bivariate and math and proposed ### ###
start_time <- Sys.time()
# model with z as input
jags.m <- jags.model(textConnection(jags_biv), 
                     data = list(z = cbind(log(thedata[,"hgt"]),
                                           log(thedata[,"wgt"])), 
                                 x1 = thedata[,"city"],
                                 x2 = thedata[,"age"],
                                 N = dim(thedata)[1]))

mAsamples <- coda.samples(jags.m, 
c("alpha0", "alpha1", "alpha2", "alpha3", "alpha4",
"gamma0", "gamma1", "gamma2", "gamma3", "gamma4", "sigma", "rho"),
n.iter = nMCMC, n.burnin=1000)
dim(mAsamples[[1]])[1]
end_time <- Sys.time()
biv_time <- end_time-start_time

#library(devtools)
#install_github("psolymos/pbapply")
library(pbapply)

# with math:
start_time <- Sys.time()
theta_samplesA <- apply(mAsamples[[1]],1, function(q) {
q["gamma1"] - 2*q["alpha1"] +
    (q["gamma3"]-2*q["alpha3"])*mean(thedata[,"age"])})  
math_est <- round(quantile(theta_samplesA, c(0.5,0.025,0.975)),3)
end_time <- Sys.time()
math_time <- (end_time-start_time) + biv_time 


# with gcomp:
start_time <- Sys.time()
x1 = thedata[,"city"]; x2 = thedata[,"age"]; S <- 2000
theta_samples <- pbapply(mAsamples[[1]], 1, create_theta_gcomp)
gcomp_est <- round(quantile(theta_samples, c(0.5,0.025,0.975)),3)
end_time <- Sys.time()
gcomp_time <- (end_time-start_time) + biv_time 





#########################
# Passive imputation strategy:



boys[,c("agecity")] <-(boys[,c("city")]*boys[,c("age")])
boys[,c("loghgt")] <- log(boys[,c("hgt")])
boys[,c("logwgt")] <- log(boys[,c("wgt")])
boys[,c("age_squared")] <- (boys[,c("age")])^2

thedata <- (boys[,c("logbmi", "loghgt", 
    "logwgt", "city", "age", "age_squared")])
dim(thedata)
start_time <- Sys.time()
dat <- thedata
init = mice(dat, maxit=0) 
head(dat)
meth = init$method
predM = init$predictorMatrix
meth
meth[c("loghgt")]="norm" 
meth[c("logwgt")]="norm" 
meth[c("age")]="norm" 
meth[c("city")]="pmm" 
meth[c("logbmi")]="~I(log( exp(logwgt)  /(( exp(loghgt)  /100)^2)))"
meth[c("age_squared")]="~I(age^2)"
meth[c("agecity")]=="~I(age*city)"
meth
n_imputations <- 50
imputed = mice(dat, method=meth, predictorMatrix=predM, m=n_imputations)
imputed_dat <- list()
for(j in 1:n_imputations){imputed_dat[[j]] <- complete(imputed,action=j)}

coda_samples <- NULL
for(j in 1:n_imputations){
  
  jags.m1 <- jags.model(textConnection(jags_univ), 
                        data = list(y = imputed_dat[[j]][,"logbmi"], 
                                    x1 = imputed_dat[[j]][,"city"],
                                    x2 = imputed_dat[[j]][,"age"],
                                    N = dim(imputed_dat[[j]])[1]))
                                                                        
  
  
  
  mAsamples <-(coda.samples(jags.m1, c("beta1", "beta3") , 
        n.iter = nMCMC, n.burnin=1000))
theta_samples <- apply(mAsamples[[1]],1, function(q) {
        q["beta1"] + q["beta3"]*mean(imputed_dat[[j]][,"age"])})  

  
  coda_samples <- c(coda_samples, theta_samples)
}
passive_samples <- coda_samples
c(mean(passive_samples), sd(passive_samples))
passive_est <- round(quantile(passive_samples, c(0.500,0.025,0.975)),3)
end_time <- Sys.time()
passive_time <- end_time-start_time


univariate_complete_est
as.numeric(univariate_complete_time[1])/60

passive_est
as.numeric(passive_time[1])/60

math_est
as.numeric(math_time[1])/60

gcomp_est
as.numeric(gcomp_time[1])/60




########################################
## Stan code for ZIKV example
########################################

Bernoulli_model <- 
  "data {
  int<lower=0> N;     // Number of observations
  int<lower=0, upper=1> y[N];   // Binary outcome data
  real<lower=0> a;               // Beta distribution shape parameter
  real<lower=0> b;               // Beta distribution shape parameter
}

parameters {
  real<lower=0, upper=1> theta;  // Probability of success
}

model {
  // Likelihood
  for (i in 1:N) {
    y[i] ~ bernoulli(theta);
  }
  
  // Priors
  theta ~ beta(a, b);
}"



BsNmN_model <- "data {
  int<lower=0> N;          
  int<lower=0, upper=N> N_z1obs;
  int<lower=0, upper=N> N_z2obs;
  int<lower=0, upper=N> N_z3obs;
  int<lower=0, upper=N> N_z1mis;
  int<lower=0, upper=N> N_z2mis;
  int<lower=0, upper=N> N_z3mis;
  vector<lower=0, upper=1>[N_z1obs] z1_obs;  // sex variable
  vector[N_z2obs] z2_obs;             // gestational age variable
  vector[N_z3obs]  z3_obs;             // head circumf. variable
  int<lower=0, upper=1> z1mis_ind[N];
  int<lower=0, upper=N> ii_z1_obs[N_z1obs];
  int<lower=0, upper=N> ii_z2_obs[N_z2obs];
  int<lower=0, upper=N> ii_z3_obs[N_z3obs];
  int<lower=0, upper=N> ii_z1_mis[N_z1mis];
  int<lower=0, upper=N> ii_z2_mis[N_z2mis];
  int<lower=0, upper=N> ii_z3_mis[N_z3mis];
}

parameters {
  vector<lower=0, upper=1>[N_z1mis] z1_mis; 
  vector[N_z2mis]  z2_mis;
  vector[N_z3mis]  z3_mis;
  real<upper=-1> kappa;
  real beta01;
  vector<lower=0>[2] zeta;    
  real beta1;
  real beta2;
  real beta3;
  simplex[2] mixweight;
  
  real mu; // mean of X
  real<lower=0> sigma; // SD of X
  real omega; // shape of X
}

transformed parameters {
  vector[N] z1;  // sex variable
  vector[N] z2;  // gestational age variable
  vector[N] z3;  // head cir variable  
  z1[ii_z1_obs] = z1_obs;
  z1[ii_z1_mis] = z1_mis;
  z2[ii_z2_obs] = z2_obs;
  z2[ii_z2_mis] = z2_mis;
  z3[ii_z3_obs] = z3_obs;
  z3[ii_z3_mis] = z3_mis;  
  vector[2] beta0;
  real loc_x; // location of X
  real gm_x; // intermediate calculation for location and scale

  gm_x = sqrt(2/pi())*omega/sqrt(1+omega^2);
  loc_x = mu - sigma*gm_x/sqrt(1-gm_x^2);
  beta0[1] = beta01;
  beta0[2] = kappa;
}

model {
  // Prior distributions for GA
  mu ~ normal(0, 0.1);
  sigma ~ inv_gamma(2, 2); 
  omega ~ normal(0, 2);
   
  // Prior distributions for HC:  
  beta0[1] ~ normal(0, 0.1);
  beta0[2] ~ normal(-2, 2)T[,-1];
  zeta[1] ~ inv_gamma(2, 2); 
  zeta[2] ~ inv_gamma(2, 2); 
  beta1 ~ normal(-0.450, 0.1);
  beta2 ~ normal(0.399, 0.1);    
  beta3 ~ normal(-0.016, 0.1); 
      
  // Likelihood
 z1_mis ~ uniform(0,1);

 z2_obs ~ skew_normal(loc_x, sigma, omega);
 z2_mis ~ skew_normal(loc_x, sigma, omega);
 
 vector[2] log_mixweight = log(mixweight);
 for (n in 1:N) {
 	vector[2] lps = log_mixweight;
// for unknown sex
   if (z1mis_ind[n]==1 ) {
		    for (k in 1:2) {
 				lps[k] += log_mix( z1[n],
                    normal_lpdf( z3[n] | 33.912 + beta0[k] + beta1 + beta2*z2[n] + 
                    beta3*pow(z2[n],2), zeta[k]),
                    normal_lpdf( z3[n] | 33.912 + beta0[k]  + beta2*z2[n]+ 
                    beta3*pow(z2[n],2), zeta[k]));
    		}	    		
  	}
// for known sex  	
	else {
			for (k in 1:2) {
		      lps[k] += normal_lpdf(z3[n] | 33.912 + beta0[k] + 
		      	beta1*z1[n] + beta2*z2[n]+ 
                    beta3*pow(z2[n],2), zeta[k]);
    		}  	 	
  	}
  	target += log_sum_exp(lps);
  }
}
"
