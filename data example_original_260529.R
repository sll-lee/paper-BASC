rm(list=ls())
#install.packages("tidysynth")
#install.packages("tmvtnorm")
library(tidysynth)
library(foreign)
library(dplyr)
library(tidyr)
library(tmvtnorm)
library(MASS)
library(MCMCpack)
library(fields)

#데이터 정리
data=read.dta("/Users/seul/Desktop/2024/BSC/realdata/repgermany.dta")
head(data)

#그외 필요한 부분 정의
vt=data %>% distinct(year) %>%pull(year)%>% as.vector() -1959
Dt=ifelse(vt-31<0, 0, 1) #intervention이 발생한 1990부터 1로 할당

#raw.y=data %>% filter(country=="West Germany") %>% pull(gdp)
data.y <- data %>% filter(country == "West Germany") %>% arrange(year)
raw.y   <- data.y$gdp
year.y <- data.y$year
pre_idx_y <- (year.y < 1990)                   # West Germany: pre-intervention 인덱스
#표준화에 사용
mean_y_pre <- mean(raw.y[pre_idx_y])
sd_y_pre   <- sd(raw.y[pre_idx_y])


table(data$year)
base.x=data %>% filter(country!="West Germany") %>%
       dplyr::select(country, year, gdp) %>% 
       pivot_wider(names_from = year, values_from = gdp)

raw.x = base.x %>% 
    dplyr::select(-country) %>%
    as.matrix()


# y <- (raw.y - mean_y_pre) / sd_y_pre

#control unit에 대해 (최종 차원은 44x16)
pre_idx=(vt<31)
t.raw.x=t(raw.x)
# mean_X_pre <- apply(t.raw.x[pre_idx, , drop=FALSE], 2, mean)
# sd_X_pre   <- apply(t.raw.x[pre_idx, , drop=FALSE], 2, sd)
# 
# x <- matrix(NA, nrow = nrow(t.raw.x), ncol = ncol(t.raw.x), 
#                 dimnames = dimnames(t.raw.x))
# 
# for (k in seq_len(ncol(x))) {
#   x[, k] <- (t.raw.x[, k] - mean_X_pre[k]) / sd_y_pre #sd_X_pre[k]
# }


#표준화 무용지물로
y= raw.y
x=t.raw.x

##chol안정화
safe_chol <- function(Sigma, base_jitter = 1e-8, max_try = 10) {
  Sigma <- (Sigma + t(Sigma)) / 2
  n <- nrow(Sigma)
  
  scale_S <- max(diag(Sigma), na.rm = TRUE)
  if (!is.finite(scale_S) || scale_S <= 0) scale_S <- 1
  
  for (k in 0:max_try) {
    jitter <- base_jitter * scale_S * 10^k
    R <- try(chol(Sigma + jitter * diag(n)), silent = TRUE)
    
    if (!inherits(R, "try-error")) {
      return(list(R = R, jitter = jitter))
    }
  }
  
  eig_min <- min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values)
  jitter <- abs(eig_min) + base_jitter * scale_S * 10^max_try
  R <- chol(Sigma + jitter * diag(n))
  
  list(R = R, jitter = jitter)
}


loglik_gamma <- function(beta, y, x, Dt, alpha.v, f, sigsq) {
  res <- as.numeric(y - x %*% beta - Dt * alpha.v - f)
  -0.5 * length(y) * log(2 * pi * sigsq) -
    0.5 * sum(res^2) / sigsq
}

#####setting
N=500000;nburn=500000
j=ncol(x)
t=length(vt)
q=1

#####function
###4) generate l
#compute posterior kernel for Metropolis
# post.dist.l = function(gp.l, dataList) {
#   tau=dataList$tau
#   vtt=dataList$vt
#   ff=dataList$f
#   alpha_l=dataList$alpha_l
#   beta_l=dataList$beta_l
#   gp.Cov=tau*exp(-0.5*rdist(vtt)^2 / gp.l)
#   # gp.Cov <- tau * exp(-0.5 * rdist(vtt)^2 / gp.l)
#   # gp.Cov <- (gp.Cov + t(gp.Cov)) / 2
#   #
#   # ch <- safe_chol(gp.Cov)
#   #
#   # gp.Cov_chol <- ch$R
#   # inv.gp.Cov <- chol2inv(gp.Cov_chol)
# 
#   if (tau < 0.001) {
#     diag(gp.Cov) = tau + 0.001
#   } else { diag(gp.Cov) = tau}
# 
#   #inverse 계산
#   gp.Cov <- (gp.Cov + t(gp.Cov))/2
#   tryCatch(
#     {
#       gp.Cov_chol <- chol(gp.Cov)
#     },
#     error = function(cond) {
#       gp.Cov <- gp.Cov + (-min(eigen(gp.Cov)$values) + 0.001)*diag(ncol(gp.Cov))
#       gp.Cov_chol <<- chol(gp.Cov)
#     }
#   )
#   inv.gp.Cov <- chol2inv(gp.Cov_chol)
# 
#   post.l=(gp.l)^(-alpha_l-1)*exp(-beta_l/gp.l)*sqrt(abs(det(gp.Cov)))^(-1)*(expm1(-0.5*t(ff)%*%inv.gp.Cov%*%ff)+1)
#   return(post.l)
# }

log_post_l <- function(gp.l, dataList) {
  tau <- dataList$tau
  vtt <- dataList$vt
  ff <- as.numeric(dataList$f)
  alpha_l <- dataList$alpha_l
  beta_l <- dataList$beta_l
  
  if (!is.finite(gp.l) || gp.l <= 0) return(-Inf)
  if (!is.finite(tau) || tau <= 0) return(-Inf)
  if (any(!is.finite(ff))) return(-Inf)
  
  gp.Cov <- tau * exp(-0.5 * rdist(vtt)^2 / gp.l)
  gp.Cov <- (gp.Cov + t(gp.Cov)) / 2
  
  ch <- try(safe_chol(gp.Cov), silent = TRUE)
  if (inherits(ch, "try-error")) return(-Inf)
  
  R <- ch$R
  
  # t(R) %*% R = gp.Cov
  z <- backsolve(R, ff, transpose = TRUE)
  quad <- sum(z^2)
  logdet <- 2 * sum(log(diag(R)))
  
  log_prior <- -(alpha_l + 1) * log(gp.l) - beta_l / gp.l
  log_gp <- -0.5 * logdet - 0.5 * quad
  
  log_prior + log_gp
}


rw_mp.l <- function(nsim.l, nburn.l, delta, dataList, initsList) {
  
  # initial value
  log.curr <- log(initsList$gp.l)
  
  # samples
  para.samples <- rep(NA_real_, nsim.l)
  acc.count <- 0
  
  for (iter.l in 1:(nsim.l + nburn.l)) {
    
    # current and proposed values
    curr <- exp(log.curr)
    log.prop <- log.curr + rnorm(1, 0, delta)
    prop <- exp(log.prop)
    
    # log acceptance ratio
    log_acc <- log_post_l(prop, dataList) -
      log_post_l(curr, dataList) +
      log.prop - log.curr
    
    # MH update
    if (is.finite(log_acc) && log(runif(1)) < log_acc) {
      log.curr <- log.prop
      
      if (iter.l > nburn.l) {
        acc.count <- acc.count + 1
      }
    }
    
    # save after burn-in
    if (iter.l > nburn.l) {
      para.samples[iter.l - nburn.l] <- exp(log.curr)
    }
    
    # optional diagnostic
    if (iter.l <= 20) {
      cat(
        "iter.l =", iter.l,
        "curr =", curr,
        "prop =", prop,
        "log_acc =", log_acc,
        "accepted =", ifelse(log.curr == log.prop, TRUE, FALSE),
        "\n"
      )
    }
  }
  
  attr(para.samples, "accept_rate") <- acc.count / nsim.l
  
  return(para.samples)
}

###4) generate u
log_post_u <- function(u, p, dataList) {
  gamma.v <- dataList$gamma.v
  alpha_u <- dataList$alpha_u
  Dt <- dataList$Dt
  alpha.v <- dataList$alpha.v
  f <- dataList$f
  x <- dataList$x
  y <- dataList$y
  sigsq <- dataList$sigsq
  
  u.j <- u[p]
  denom <- sum(gamma.v * u)
  
  if (!is.finite(u.j) || u.j <= 0) return(-Inf)
  if (!is.finite(denom) || denom <= 0) return(-Inf)
  if (!is.finite(sigsq) || sigsq <= 0) return(-Inf)
  
  beta.tem <- gamma.v * u / denom
  
  res <- as.numeric(y - x %*% beta.tem - Dt * alpha.v - f)
  
  log_prior <- (alpha_u - 1) * log(u.j) - alpha_u * u.j
  log_lik <- -0.5 * sum(res^2) / sigsq
  
  log_prior + log_lik
}
# 
# post.dist.u = function(u, p, dataList) {
#   gamma.v = dataList$gamma.v
#   alpha_u=dataList$alpha_u
#   Dt = dataList$Dt
#   alpha.v = dataList$alpha.v
#   f = dataList$f
#   x = dataList$x
#   y = dataList$y
#   sigsq = dataList$sigsq
#   u.tem=u
#   u.j=u.tem[p]
#   
#   beta.tem=rep(0,j)
#   for (k in 1:j) {beta.tem[k]=gamma.v[k]*u.tem[k] / sum(gamma.v*u.tem)}
#   post.u=u.j^(alpha_u-1)*exp(-alpha_u*u.j)*exp(-sum((y-x%*%beta.tem-Dt*alpha.v-f)^2)/(2*sigsq))
#   return(post.u)
# }

# rw_mp.u=function(nsim.u, nburn.u, delta, p, dataList, initsList){
#   #initial value
#   log.curr=log(initsList$u[p])
#   
#   #start iterations
#   para.samples=rep(0,nsim.u)
#   for (iter.u in 1:(nsim.u+nburn.u)){
#     z=rnorm(1,0,1)
#     curr=exp(log.curr)
#     log.prop=z*delta+log.curr
#     prop=exp(log.prop)
#     u.curr=initsList$u
#     u.curr[p]=curr
#     u.prop=initsList$u
#     u.prop[p]=prop
#     
#     select.alpha.u=log(post.dist.u(u.prop, p, dataList))-log(post.dist.u(u.curr, p, dataList)) +log(prop)-log(curr)
#     
#     if(is.finite(select.alpha.u)){
#       if(log(runif(1))<select.alpha.u) {para.next = log.prop} else {para.next = log.curr}
#     } else {para.next=log.curr}
#     
#     log.curr=para.next
#     if (iter.u>nburn.u) para.samples[iter.u-nburn.u]=exp(para.next)
#     
#   }
#   
#   return(para.samples)
# }


german.chain = function(seed) {

## start gibbs sampling
f.sample = matrix(0, N,t)
sig.sample = rep(0, N)
tau.sample = rep(0, N)
l.sample = rep(0, N)
w.sample =rep(0, N)
gamma.sample = matrix(0, N, j)
gamma.post=rep(0,j)
u.sample = matrix(0, N, j)
a.sample = matrix(0, N, q)

#initial point
#prior parameter
alpha_sig <- 10
beta_sig  <- 5000

alpha_tau <- 3
beta_tau <- 2*10000

alpha_l=3
beta_l=1000

alpha_u = 2.5

# initial point :240829 수정
set.seed(seed)
tau = rinvgamma(n=1, alpha_tau, beta_tau)
#tau = (beta_tau) / (alpha_tau - 1) 
print(tau)
lscale = rinvgamma(n=1, alpha_l, beta_l)
#lscale = beta_l / (alpha_l - 1) 
print(lscale)
Cov.l = matrix(0, nrow=t, ncol=t)
for (i in 1:t){
  for (k in 1:t){
    Cov.l[i,k]=exp(-0.5*(vt[i]-vt[k])^2 / lscale)
  }
}
mu=rep(0,t)
f = mvrnorm(n=1, mu, tau*Cov.l)
print(f)
w = runif(1,0,1)

#gamma.v <- rbinom(j, 1, w)
gamma.v = rep(1,j)
if (sum(gamma.v) == 0) {
  gamma.v[sample(1:j, 1)] <- 1
}
u = rgamma(j, alpha_u, scale=1/alpha_u)
print(u)
beta=rep(0,j)
for (k in 1:j) {beta[k]=gamma.v[k]*u[k] / sum(gamma.v*u)}
sigsq=rinvgamma(1, alpha_sig, beta_sig)
print(sigsq)


sig.a = rinvgamma(1, 3, 2*1000^2)
mu_a = rep(-0,q)
#alpha.v = mvrnorm(n=1, mu_a, sig.a*diag(q))
alpha.v = rtmvnorm(n = 1, mean = mu_a, sigma = sig.a * diag(q), upper = rep(0, q)) #truncated 240923



# Initialize delta.v for u
#delta.v <- ifelse(gamma.v == 1, 1, 2)
delta.v=rep(3,j)
delta.fixed <- rep(FALSE, j)
# Adaptive Metropolis-Hastings for u
adapt.iter <- 1000
accept.count <- rep(0, j)
accept.count.u <- rep(0, j)
proposal.count.u <- rep(0, j)
target_acceptance_rate <- 0.23
adapt_rate <- 0.05
u_update_freq = 5

for (iter in 1:(N+nburn)){
  
  ## 1) generate f
  
  # GP covariance
  Cov.f <- tau * exp(-0.5 * rdist(vt)^2 / lscale)
  Cov.f <- (Cov.f + t(Cov.f)) / 2
  
  # Cholesky of Cov.f
  ch_Cov <- safe_chol(Cov.f)
  Cov.f_chol <- ch_Cov$R
  
  # inverse covariance
  inv.Cov.f <- chol2inv(Cov.f_chol)
  inv.Cov.f <- (inv.Cov.f + t(inv.Cov.f)) / 2
  
  # beta
  if (sum(gamma.v * u) <= 0) {
    stop(paste("Empty donor set in f update at iter", iter))
  }
  
  beta <- gamma.v * u / sum(gamma.v * u)
  
  # posterior precision of f
  Q.f <- inv.Cov.f + diag(1 / sigsq, t)
  Q.f <- (Q.f + t(Q.f)) / 2
  
  # Cholesky of posterior precision
  ch_Q <- safe_chol(Q.f)
  R_Q <- ch_Q$R
  
  # posterior mean: mu.f = solve(Q.f, (y - X beta - D alpha) / sigsq)
  b.f <- as.numeric(y - x %*% beta - Dt * alpha.v) / sigsq
  
  tmp <- forwardsolve(t(R_Q), b.f)
  mu.f <- backsolve(R_Q, tmp)
  
  # sample from N(mu.f, Q.f^{-1})
  z <- rnorm(t)
  f <- as.vector(mu.f + backsolve(R_Q, z))
  
# 
#   #1) generate f
# 
#   Cov.f = tau*exp(-0.5*rdist(vt)^2 / lscale)
#   Cov.f <- (Cov.f + t(Cov.f))/2
#   tryCatch(
#     {
#       Cov.f_chol <- chol(Cov.f)
#     },
#     error = function(cond) {
#       Cov.f <- Cov.f + (-min(eigen(Cov.f)$values) + 0.001)*diag(ncol(Cov.f))
#       Cov.f_chol <<- chol(Cov.f)
#     }
#   )
#   inv.Cov.f <- chol2inv(Cov.f_chol)
# 
#   V.inv = inv.Cov.f+sigsq^(-1)*diag(1,t)
# 
# 
#   #V.f 재정의:240815
#   V.inv <- (V.inv + t(V.inv))/2
#   tryCatch(
#     {
#       Cov.f_chol <- chol(V.inv)
#     },
#     error = function(cond) {
#       V.inv <- V.inv + (-min(eigen(V.inv)$values) + 0.001)*diag(ncol(V.inv))
#       Cov.f_chol <<- chol(V.inv)
#     }
#   )
#   V.f <- chol2inv(Cov.f_chol)
# 
# 
#   #V.f=solve(V.inv):240815 비활성화
# 
#   beta=rep(0,j)
#   if (sum(gamma.v*u) != 0) {
#     for (k in 1:j) {beta[k]=gamma.v[k]*u[k] / sum(gamma.v*u)}
#   } else {
#     beta=rep(0,j)
#   }
#   mu.f = V.f%*%(y-x%*%beta-Dt*alpha.v)/sigsq
#   f = as.vector(mvrnorm(1, mu.f, V.f))

  
  #2) generate sigma^2
  beta=rep(0,j)
  
  if (sum(gamma.v*u) != 0) {
    for (k in 1:j) {beta[k]=gamma.v[k]*u[k] / sum(gamma.v*u)}
  } else {
    beta=rep(0,j)
  }
  
  alpha_sig.post=alpha_sig+t/2
  beta_sig.post=beta_sig+1/2*sum((y-x%*%beta-Dt*alpha.v-f)^2)
  sigsq=rinvgamma(1, alpha_sig.post, beta_sig.post)
  
  
  #3) generate tau^2
  alpha_tau.post=alpha_tau+t/2
  Cov.l = exp(-0.5*rdist(vt)^2 / lscale)
  
  Cov.l <- (Cov.l + t(Cov.l))/2
  tryCatch(
    {
      Cov.l_chol <- chol(Cov.l)
    },
    error = function(cond) {
      Cov.l <- Cov.l + (-min(eigen(Cov.l)$values) + 0.001)*diag(ncol(Cov.l))
      Cov.l_chol <<- chol(Cov.l)
    }
  )
  inv.Cov.l <- chol2inv(Cov.l_chol)
  
  beta_tau.post=beta_tau+0.5*t(f)%*%inv.Cov.l%*%f
  tau=rinvgamma(1, alpha_tau.post, beta_tau.post)
  
  #4) generate l
  #prior & data
  dataList.l=list(tau=tau, vt=vt, f=f, alpha_l=alpha_l, beta_l=beta_l)
  delta=0.5
  initsList=list(gp.l=lscale)
  lscale=rw_mp.l(1,0, delta,dataList.l,initsList)
  
  #new:250826
#  if (lscale < 1)  lscale <- 1
#  if (lscale > 64) lscale <- 64
  
  #5) generate w
  alpha_w.post=sum(gamma.v)+1
  beta_w.post=sum(rep(1,j)-gamma.v)+1
  w=rbeta(1,alpha_w.post,beta_w.post)
  
  
  #6) generate gamma_j
  #  while (TRUE){
  
  loglik_gamma <- function(beta, y, x, Dt, alpha.v, f, sigsq) {
    res <- as.numeric(y - x %*% beta - Dt * alpha.v - f)
    -0.5 * length(y) * log(2 * pi * sigsq) -
      0.5 * sum(res^2) / sigsq
  }
  
  for (p in 1:j) {
    
    ## 6-1) gamma_p = 1
    gamma.tem1 <- gamma.v
    gamma.tem1[p] <- 1
    
    denom1 <- sum(gamma.tem1 * u)
    
    if (denom1 <= 0 || !is.finite(denom1)) {
      log_rj_1 <- -Inf
    } else {
      beta.tem1 <- gamma.tem1 * u / denom1
      log_rj_1 <- log(w) + loglik_gamma(
        beta = beta.tem1,
        y = y,
        x = x,
        Dt = Dt,
        alpha.v = alpha.v,
        f = f,
        sigsq = sigsq
      )
    }
    
    ## 6-2) gamma_p = 0
    gamma.tem0 <- gamma.v
    gamma.tem0[p] <- 0
    
    denom0 <- sum(gamma.tem0 * u)
    
    if (denom0 <= 0 || !is.finite(denom0)) {
      # empty donor set은 허용하지 않음
      log_rj_0 <- -Inf
    } else {
      beta.tem0 <- gamma.tem0 * u / denom0
      log_rj_0 <- log(1 - w) + loglik_gamma(
        beta = beta.tem0,
        y = y,
        x = x,
        Dt = Dt,
        alpha.v = alpha.v,
        f = f,
        sigsq = sigsq
      )
    }
    
    ## stable Bernoulli probability
    bin.p <- plogis(log_rj_1 - log_rj_0)
    
    if (!is.finite(bin.p) || bin.p < 0 || bin.p > 1) {
      print(list(
        p = p,
        log_rj_1 = log_rj_1,
        log_rj_0 = log_rj_0,
        sigsq = sigsq,
        sum_gamma_current = sum(gamma.v),
        denom1 = denom1,
        denom0 = denom0
      ))
      stop("Bad gamma probability")
    }
    
    gamma.v[p] <- rbinom(1, 1, bin.p)
  }
  # for (p in 1:j){
  #   gamma.tem1=gamma.v
  #   #6-1) r_j=1
  #   gamma.tem1[p]=1
  #   beta.tem1=rep(0,j)
  #   for (k in 1:j) {beta.tem1[k]=gamma.tem1[k]*u[k] / sum(gamma.tem1*u)}
  #   value1=(2*pi*sigsq)^(-t/2)*exp(-sum((y-x%*%beta.tem1-Dt*alpha.v-f)^2)/(2*sigsq))
  #   rj_1=w*value1
  #   
  #   #6-2)r_j=0
  #   gamma.tem0=gamma.v
  #   gamma.tem0[p]=0
  #   beta.tem0=rep(0,j)
  #   if (sum(gamma.tem0*u) != 0) {
  #     for (k in 1:j) {beta.tem0[k]=gamma.tem0[k]*u[k] / sum(gamma.tem0*u)}
  #   } else {
  #     beta.tem0=rep(0,j)
  #   }
  #   value0=(2*pi*sigsq)^(-t/2)*exp(-sum((y-x%*%beta.tem0-Dt*alpha.v-f)^2)/(2*sigsq))
  #   rj_0=(1-w)*value0
  #   
  #   bin.p=(rj_1)/(rj_1+rj_0)
  #   gamma.v[p]=rbinom(1,1,bin.p)
  # }
  
  
  # #7) generate u  
  # #prior & data
  # dataList.u=list(gamma.v=gamma.v, alpha_u=alpha_u, Dt=Dt, alpha.v=alpha.v,f=f, x=x, y=y, sigsq=sigsq)
  # initsList=list(u=u)
  # 
  # ##iter 별도 추가 : 240806
  # 
  # for (freq in 1:u_update_freq) {
  #   for (p in 1:j) {
  #     delta <- delta.v[p]
  #     u.curr = initsList$u[p]
  #     u.prop <- initsList$u
  #     u.prop[p] <-exp(log(u.prop[p])+rnorm(1,0,1)*delta)
  #     
  #     log_accept_prob <- log(post.dist.u(u.prop, p, dataList.u)) - log(post.dist.u(initsList$u, p, dataList.u)) +log(u.prop[p])-log(u.curr)
  #     
  #     if (is.finite(log_accept_prob) && log(runif(1)) < log_accept_prob) {
  #       initsList$u <- u.prop
  #       #accept.count[p] <- accept.count[p] + 1
  #     } else {
  #       # Do nothing, retain the current value of initsList$u
  #     }
  #   }
  # }
  # 
  # 
  # u=initsList$u
  
  
  #7) generate u
  dataList.u <- list(
    gamma.v = gamma.v,
    alpha_u = alpha_u,
    Dt = Dt,
    alpha.v = alpha.v,
    f = f,
    x = x,
    y = y,
    sigsq = sigsq
  )
  
  initsList <- list(u = u)
  
  for (freq in 1:u_update_freq) {
    for (p in 1:j) {
      
      delta <- delta.v[p]
      
      u.curr <- initsList$u[p]
      u.prop <- initsList$u
      u.prop[p] <- exp(log(u.prop[p]) + rnorm(1, 0, 1) * delta)
      
      log_accept_prob <-
        log_post_u(u.prop, p, dataList.u) -
        log_post_u(initsList$u, p, dataList.u) +
        log(u.prop[p]) - log(u.curr)
      
      # proposal count
      if (iter <= nburn) {
        proposal.count.u[p] <- proposal.count.u[p] + 1
      }
      
      if (is.finite(log_accept_prob) && log(runif(1)) < log_accept_prob) {
        initsList$u <- u.prop
        
        # acceptance count
        if (iter <= nburn) {
          accept.count.u[p] <- accept.count.u[p] + 1
        }
      }
    }
  }
  
  u <- initsList$u
  
  #8) generate a
  W.inv = sig.a^(-1)*diag(1,q)+sigsq^(-1)*t(Dt)%*%Dt
  W.a = solve(W.inv)
  
  W.inv <- (W.inv + t(W.inv))/2
  tryCatch(
    {
      W.inv_chol <- chol(W.inv)
    },
    error = function(cond) {
      W.inv <- W.inv + (-min(eigen(W.inv)$values) + 0.001)*diag(ncol(W.inv))
      W.inv_chol <<- chol(W.inv)
    }
  )
  W.a <- chol2inv(W.inv_chol)
  
  beta=rep(0,j)
  if (sum(gamma.v*u) != 0) {
    for (k in 1:j) {beta[k]=gamma.v[k]*u[k] / sum(gamma.v*u)}
  } else {
    beta=rep(0,j)
  }
  
  mu.a = c(sigsq^(-1)*W.a%*%t(Dt)%*%(y-x%*%beta-f)) #c() 추가
#  alpha.v = as.vector(mvrnorm(1, mu.a, W.a))
  alpha.v = as.vector(rtmvnorm(1, mu.a, W.a, upper = rep(0,q))) #truncated 240923

  
  ## diagnostic check
  if (!is.finite(sigsq) || sigsq <= 0) stop(paste("Bad sigsq at iter", iter))
  if (!is.finite(tau) || tau <= 0) stop(paste("Bad tau at iter", iter))
  if (any(!is.finite(f))) stop(paste("Bad f at iter", iter))
  if (any(!is.finite(u)) || any(u <= 0)) stop(paste("Bad u at iter", iter))
  if (any(!is.finite(alpha.v))) stop(paste("Bad alpha.v at iter", iter))
  if (any(!is.finite(beta))) stop(paste("Bad beta at iter", iter))
  
  if (iter %% 100 == 0) {
    res <- as.numeric(y - x %*% beta - Dt * alpha.v - f)
    cat(
      "iter =", iter,
      "sigsq =", sigsq,
      "tau =", tau,
      "rss =", sum(res^2),
      "max|res| =", max(abs(res)),
      "sum_gamma =", sum(gamma.v),
      "\n"
    )
  }
  
  if (iter>nburn){
    f.sample[iter-nburn,]=f
    sig.sample[iter-nburn]=sigsq
    tau.sample[iter-nburn]=tau
    l.sample[iter-nburn]=lscale
    w.sample[iter-nburn]=w
    gamma.sample[iter-nburn,]=gamma.v
    u.sample[iter-nburn,]=u
    a.sample[iter-nburn,]=alpha.v
  }
  
  # if (iter > nburn && (iter - nburn) %% adapt.iter == 0) {
  #   met.u <- mcmc(u.sample[1:(iter - nburn), ])  # 현재까지의 샘플로 mcmc 객체 생성
  #   acceptance.rate <- 1 - rejectionRate(met.u)
  #   print(acceptance.rate)
  #   print(delta.v)
  #   #acceptance.rate <- accept.count / (adapt.iter * u_update_freq)
  #   for (p in 1:j) {
  #     if (!delta.fixed[p]) {  # delta가 고정되지 않은 경우에만 조정
  #       delta_min <- 0.05
  #       delta_max <- 6
  #       
  #       if (acceptance.rate[p] < target_acceptance_rate * 0.95) {
  #         delta.v[p] <- max(delta.v[p] * (1 - adapt_rate), delta_min)
  #       } else if (acceptance.rate[p] > target_acceptance_rate * 1.05) {
  #         delta.v[p] <- min(delta.v[p] * (1 + adapt_rate), delta_max)
  #       } else {
  #         delta.fixed[p] <- TRUE
  #       }
  #       # 
  #       # if (acceptance.rate[p] < target_acceptance_rate * 0.95) {
  #       #   delta.v[p] <- max(delta.v[p] * (1 - adapt_rate),1)
  #       # } else if (acceptance.rate[p] > target_acceptance_rate * 1.05) {
  #       #   delta.v[p] <- min(delta.v[p] * (1 + adapt_rate),6)
  #       # } else {
  #       #   delta.fixed[p] <- TRUE  # 채택률이 목표 수준에 도달하면 delta 고정
  #       # }
  #     }
  #   }
  # }
  
  if (iter <= nburn && iter %% adapt.iter == 0) {
    
    acceptance.rate <- accept.count.u / pmax(proposal.count.u, 1)
    
    print(acceptance.rate)
    print(delta.v)
    
    for (p in 1:j) {
      
      delta_min <- 0.05
      delta_max <- 6
      
      if (acceptance.rate[p] < target_acceptance_rate * 0.95) {
        delta.v[p] <- max(delta.v[p] * (1 - adapt_rate), delta_min)
        
      } else if (acceptance.rate[p] > target_acceptance_rate * 1.05) {
        delta.v[p] <- min(delta.v[p] * (1 + adapt_rate), delta_max)
      }
    }
    
    # 다음 adapt.iter 구간을 위해 count 초기화
    accept.count.u[] <- 0
    proposal.count.u[] <- 0
  }
  
  print(iter)
  print(lscale)
}
return(list(
  y=y,
  x=x,
  f.sample = f.sample,
  sig.sample = sig.sample,
  tau.sample = tau.sample,
  l.sample = l.sample,
  w.sample = w.sample,
  gamma.sample = gamma.sample,
  u.sample = u.sample,
  a.sample = a.sample
))
}


seeds = c(103,200,305)
#seeds = c(200,305)
#seeds=c(200,300)
#seeds=200
results = list()

#results[[3]] = german.chain(305)


for (seed in seeds) {
  results[[seed%/%100]] <- german.chain(seed)
}

#별도 저장 / 다시 불러올 때는 load() 활용
save(results, file = "/Users/seul/Desktop/2024/BSC/realdata/260518 revision/B-our_0529.RData")

load(file = "/Users/seul/Desktop/2024/BSC/realdata/260518 revision/B-our_0529.RData")


met=mcmc(results[[2]]$l.sample)
print(1-rejectionRate(met))

met=mcmc(results[[2]]$u.sample)
print(1-rejectionRate(met))

plot(results[[2]]$l.sample,type="l")
mean(results[[2]]$tau.sample)
mean(results[[2]]$l.sample)

for (i in 1:3){
  met=mcmc(results[[i]]$l.sample)
  print(1-rejectionRate(met))
}

for (i in 1:3){
  met.u=mcmc(results[[i]]$u.sample)
  #met.u=mcmc(u.sample)
  print(1-rejectionRate(met.u))
}

for (i in 1:3){
  results[[i]]$beta.sample = matrix(0, N, j)
  for (beta.row in 1:N) {
    for (k in 1:j) {results[[i]]$beta.sample[beta.row,k]=results[[i]]$gamma.sample[beta.row, k]*results[[i]]$u.sample[beta.row, k] / sum(results[[i]]$gamma.sample[beta.row,]*results[[i]]$u.sample[beta.row,])}
  }
}




##################gelman


sample.names = c("f.sample","sig.sample","tau.sample","l.sample","w.sample", "u.sample", "a.sample")

# 결과 저장 리스트
gelman.summary = list()

# 각 샘플에 대해 겔만 상수 계산
for (sample.name in sample.names) {
  # 각 체인에서 해당 샘플 추출 및 mcmc 객체로 변환
  mcmc.chain1 = (results[[1]][[sample.name]])
  mcmc.chain2 = (results[[2]][[sample.name]])
  mcmc.chain3 = (results[[3]][[sample.name]])
  #  mcmc.chain4 = (results[[4]][[sample.name]])
  #  mcmc.chain5 = (results[[5]][[sample.name]])
  
  #thin개 중 하나씩 샘플링하여 새로운 체인 생성
  thin = 1
  if (is.matrix(mcmc.chain1)) {
    new.chain1 <- mcmc(mcmc.chain1[seq(1, nrow(mcmc.chain1), by=thin), ])
    new.chain2 <- mcmc(mcmc.chain2[seq(1, nrow(mcmc.chain2), by=thin), ])
    new.chain3 <- mcmc(mcmc.chain3[seq(1, nrow(mcmc.chain3), by=thin), ])
    #    new.chain4 <- mcmc(mcmc.chain4[seq(1, nrow(mcmc.chain4), by=thin), ])
    #   new.chain5 <- mcmc(mcmc.chain5[seq(1, nrow(mcmc.chain5), by=thin), ])
  } else {
    new.chain1 = mcmc(mcmc.chain1[seq(1, length(mcmc.chain1), by=thin)])
    new.chain2 = mcmc(mcmc.chain2[seq(1, length(mcmc.chain2), by=thin)])
    new.chain3 = mcmc(mcmc.chain3[seq(1, length(mcmc.chain3), by=thin)])
    #  new.chain4 = mcmc(mcmc.chain4[seq(1, length(mcmc.chain4), by=thin)])
    # new.chain5 = mcmc(mcmc.chain5[seq(1, length(mcmc.chain5), by=thin)])
  }
  
  # mcmc.list 생성
  mcmc.list = mcmc.list(new.chain1, new.chain2, new.chain3)
  #, new.chain4, new.chain5)
  
  # 겔만 상수 계산
  #print(gelman.plot(mcmc.list))
  gelman_diag = gelman.diag(mcmc.list)
  
  # 결과 저장
  gelman.summary[[sample.name]] <- gelman_diag
}

gelman.summary

beta.gelman <- list()

for (k in 1:j) {
  beta.gelman[[k]] <- gelman.diag(mcmc.list(
    mcmc(results[[1]]$beta.sample[, k]),
    mcmc(results[[2]]$beta.sample[, k]),
    mcmc(results[[3]]$beta.sample[, k])
  ))
  
  cat("\nBeta", k, "\n")
  print(beta.gelman[[k]])
}


######################predict y


u.sample=results[[2]]$u.sample
sig.sample=results[[2]]$sig.sample
beta.sample=results[[2]]$beta.sample
a.sample=results[[2]]$a.sample
f.sample=results[[2]]$f.sample
y=results[[2]]$y
x=results[[2]]$x

eps.sample = matrix(0, N, t)
for (i in 1:N){
  eps.sample[i,] = rnorm(t, mean=0, sd=sig.sample[i])
}

y.sample = matrix(0, N, t)
for (i in 1:N){
  y.sample[i,]=x%*%beta.sample[i,]+Dt*0 +f.sample[i,]
  #+eps.sample[i]
}

y.ci = apply(y.sample, 2, quantile, probs = c(0.025, 0.5, 0.975))

y.predict = data.frame(
  vt=vt,
  y_mean = y.ci[2,],
  ymin = y.ci[1,],
  ymax = y.ci[3,]
)

y.true = data.frame(vt=vt, y = y)

library(ggplot2)

ggplot(y.predict, aes(x = vt, y = y_mean)) +
  # ribbon: ymin/zmin, ymax/zmax
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "skyblue", alpha = 0.3) +
  # z_mean 선 그래프
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 2) +
  scale_y_continuous(limits = c(0, 35000))



#Table5 : B-our
round(colMeans(beta.sample),2)
round(apply(beta.sample, 2, sd),2)

#Table5 : B-MV
round(colMeans(beta.sample2),2)
round(apply(beta.sample2, 2, sd),2)


save(y.predict, file = "/Users/seul/Desktop/2024/BSC/realdata/260518 revision/comparison/260529_realdata_B_our.RData")

#load("/Users/seul/Desktop/2024/BSC/realdata/realdata_standard_included_f_250615.RData")
write.csv(y.predict, file = "/Users/seul/Desktop/2024/BSC/realdata/260518 revision/comparison/z.predict.95_0529.csv", row.names = FALSE)
#write.csv(z.predict2, file = "/Users/seul/Desktop/2024/BSC/realdata/z.predict2.95_250630.csv", row.names = FALSE)


#--------------------------------

library(ggplot2)
#plot1
ggplot() +
#  geom_ribbon(data = y.predict, aes(x=vt, ymin = ymin, ymax = ymax), fill = "blue", alpha = 0.2)+
  #  geom_ribbon(data = y.predict2, aes(x=vt, ymin = ymin, ymax = ymax), fill = "grey", alpha = 0.7)+
  geom_line(data = y.predict, aes(x = vt, y = y_mean), color = "red", linetype = "dashed")+
  #geom_line(data = y.predict, aes(x = vt, y = y_mean), color = "blue", linetype = "dashed")+
  geom_line(data = y.true, aes(x = vt, y = y), color = "black", linewidth = 0.5)+
  labs(x="time", y="y")



#####################################################################
######################250216 traceplot###############################
#####################################################################
library(ggplot2)
library(reshape2)
setwd("/Users/seul/Desktop/2024/BSC/realdata/260518 revision")
# (1) 숫자에 접미사(st, nd, rd, th)를 붙이는 함수
ordinal_suffix <- function(num) {
  suffix <- ifelse(num %% 10 == 1 & num %% 100 != 11, "st",
                   ifelse(num %% 10 == 2 & num %% 100 != 12, "nd",
                          ifelse(num %% 10 == 3 & num %% 100 != 13, "rd", "th")))
  return(paste0(num, suffix))
}

plot_all_traces <- function(results, thin = 100,
                            selected_u = c(1, 10, 2, 7),
                            selected_f = c(1, 10, 40),
                            output_dir = "plot/") {
  dir.create(output_dir, showWarnings = FALSE)  # 폴더 생성
  
  trace_list <- list()
  
  #### 1. 스칼라 파라미터들
  scalar_params <- c("sig.sample", "tau.sample", "l.sample", "w.sample", "a.sample")
  for (param in scalar_params) {
    chain_data <- lapply(seq_along(results), function(chain_id) {
      iter_idx <- seq(1, nrow(as.matrix(results[[chain_id]][[param]])), by = thin)
      df <- data.frame(iter      = iter_idx,
                       value     = results[[chain_id]][[param]][iter_idx],
                       chain     = as.factor(chain_id),
                       parameter = param)
      return(df)
    })
    df_param <- do.call(rbind, chain_data)
    trace_list[[param]] <- df_param
  }
  
  #### 2. u.sample
  if ("u.sample" %in% names(results[[1]])) {
    chain_data <- lapply(seq_along(results), function(chain_id) {
      mat <- results[[chain_id]][["u.sample"]]
      iter_idx <- seq(1, nrow(mat), by = thin)
      available_u <- ncol(mat)
      sel_u <- selected_u[selected_u <= available_u]
      df <- as.data.frame(mat[iter_idx, sel_u, drop = FALSE])
      colnames(df) <- sapply(sel_u, ordinal_suffix)  # "1st", "10th", ...
      df$iter <- iter_idx
      df$chain <- as.factor(chain_id)
      
      df_melt <- melt(df, id.vars = c("iter", "chain"),
                      variable.name = "sub_index", value.name = "value")
      # factor 레벨 지정: sub_index 순서를 selected_u 순서대로 고정
      df_melt$sub_index <- factor(
        df_melt$sub_index,
        levels = sapply(sel_u, ordinal_suffix)  # 예: c("1st", "10th", "2nd", "7th")
      )
      # parameter도 factor로 지정 (이후 최종 순서에서 재정렬할 예정)
      df_melt$parameter <- paste0("u (", df_melt$sub_index, " control unit)")
      
      return(df_melt[, c("iter", "value", "chain", "parameter")])
    })
    df_u <- do.call(rbind, chain_data)
    trace_list[["u.sample"]] <- df_u
  }
  
  #### 3. f.sample
  if ("f.sample" %in% names(results[[1]])) {
    chain_data <- lapply(seq_along(results), function(chain_id) {
      mat <- results[[chain_id]][["f.sample"]]
      iter_idx <- seq(1, nrow(mat), by = thin)
      available_f <- ncol(mat)
      sel_f <- selected_f[selected_f <= available_f]
      df <- as.data.frame(mat[iter_idx, sel_f, drop = FALSE])
      colnames(df) <- as.character(sel_f)
      df$iter <- iter_idx
      df$chain <- as.factor(chain_id)
      
      df_melt <- melt(df, id.vars = c("iter", "chain"),
                      variable.name = "sub_index", value.name = "value")
      # "f.sample (t=...)" 를 "f (t=...)" 로 저장
      df_melt$parameter <- paste0("f (t=", df_melt$sub_index, ")")
      return(df_melt[, c("iter", "value", "chain", "parameter")])
    })
    df_f <- do.call(rbind, chain_data)
    trace_list[["f.sample"]] <- df_f
  }
  
  #### 4. 합치기
  all_trace <- do.call(rbind, trace_list)
  
  # --- (A) 부분 문자열 치환 ---
  # "f.sample (t=1)" -> "f (t=1)"
  all_trace$parameter <- gsub("f.sample", "f", all_trace$parameter)
  # "u.sample (10th control unit)" -> "u (10th control unit)"
  all_trace$parameter <- gsub("u.sample", "u", all_trace$parameter)
  
  # 스칼라 파라미터 치환
  all_trace$parameter <- gsub("sig.sample", "sigma^2", all_trace$parameter)
  all_trace$parameter <- gsub("w.sample",   "eta",       all_trace$parameter)
  all_trace$parameter <- gsub("a.sample",   "alpha",     all_trace$parameter)
  all_trace$parameter <- gsub("tau.sample", "tau",       all_trace$parameter)
  all_trace$parameter <- gsub("l.sample",   "kappa",     all_trace$parameter)
  
  # --- (B) 전체 facet 순서 지정 ---
  # 원하는 최종 순서를 벡터로 지정
  final_levels <- c(
    "f (t=1)", "f (t=10)", "f (t=40)",
    "tau", "eta", "alpha", "sigma^2", "kappa",
    "u (1st control unit)", "u (10th control unit)", "u (2nd control unit)", "u (7th control unit)"
  )
  # 실제 데이터에 존재하는 레벨만 유지
  final_levels <- intersect(final_levels, unique(all_trace$parameter))
  all_trace$parameter <- factor(all_trace$parameter, levels = final_levels)
  
  # --- (C) labeller 함수 ---
  scalar_labels <- c("sigma^2", "tau", "eta", "alpha", "kappa")
  custom_labeller <- function(x) {
    x <- as.character(x)
    lapply(x, function(lbl) {
      if (lbl %in% scalar_labels) {
        parse(text = lbl)[[1]]
      } else {
        lbl
      }
    })
  }
  
  n_facets <- length(unique(all_trace$parameter))
  n_rows <- ceiling(n_facets / 3)
  total_height <- n_rows * 2
  
  # facet_wrap
  p <- ggplot(all_trace, aes(x = iter, y = value, color = chain)) +
    geom_line(alpha = 0.7) +
    scale_x_continuous(
      breaks = c(100000, 300000, 500000),
      labels = function(x) format(x, big.mark = ",", scientific = FALSE)
    )+
    facet_wrap(~ parameter, scales = "free_y", ncol = 3, labeller = as_labeller(custom_labeller)) +
    labs(title = "", x = "Iteration", y = "Value") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 12),         # 축 눈금 크기
      legend.title = element_text(size = 16),        # 범례 제목 크기
      legend.text = element_text(size = 14),         # 범례 항목 크기
      strip.text = element_text(size = 18, face = "bold")
    )
  
  
  # 저장
  dir.create(output_dir, showWarnings = FALSE)
  
  filename <- file.path(output_dir, "trace_plots_realdata_0529")
  
  # PNG
  ggsave(paste0(filename, ".png"), plot = p, width = 12, height = total_height, dpi = 150)
  
  # PDF
  ggsave(paste0(filename, ".pdf"), plot = p, width = 12, height = total_height)
  
  # EPS (권장: cairo_ps 사용 → 글꼴/투명도 문제 방지)
  ggsave(paste0(filename, ".eps"), plot = p, width = 12, height = total_height,
         device = cairo_ps, fallback_resolution = 600)
  
  return(p)
}

# 사용 예시
plot_all_traces(results, thin = 100)


#####################################################################
#################250216 marginal density#############################
#####################################################################
library(ggplot2)
library(reshape2)

# 숫자를 받아서 "1st", "2nd", "3rd", ... 로 변환하는 함수
ordinal_suffix <- function(num) {
  suffix <- ifelse(num %% 10 == 1 & num %% 100 != 11, "st",
                   ifelse(num %% 10 == 2 & num %% 100 != 12, "nd",
                          ifelse(num %% 10 == 3 & num %% 100 != 13, "rd", "th")))
  return(paste0(num, suffix))
}

plot_all_densities <- function(results, thin = 100, 
                               selected_u = c(1, 10, 2, 7),    
                               selected_f = c(1, 10, 25, 50),  
                               output_dir = "plot/") {
  dir.create(output_dir, showWarnings = FALSE)  # 폴더 생성
  
  density_list <- list()
  
  #### 1. 스칼라 파라미터들: sig.sample, w.sample, a.sample
  scalar_params <- c("sig.sample", "tau.sample","l.sample","w.sample", "a.sample")
  for (param in scalar_params) {
    chain_data <- lapply(seq_along(results), function(chain_id) {
      param_vals <- results[[chain_id]][[param]]
      thin_idx <- seq(1, length(param_vals), by = thin)
      df <- data.frame(value = as.vector(param_vals[thin_idx]),
                       chain = as.factor(chain_id),
                       parameter = param)
      return(df)
    })
    df_param <- do.call(rbind, chain_data)
    density_list[[param]] <- df_param
  }
  
  #### 2. u.sample: 선택한 열만 추출 (모든 iteration에서 thin 적용)
  if ("u.sample" %in% names(results[[1]])) {
    chain_data <- lapply(seq_along(results), function(chain_id) {
      mat <- results[[chain_id]][["u.sample"]]
      available_u <- ncol(mat)
      sel_u <- selected_u[selected_u <= available_u]
      thin_idx <- seq(1, nrow(mat), by = thin)
      
      df_u <- do.call(rbind, lapply(seq_along(sel_u), function(i) {
        data.frame(value = mat[thin_idx, sel_u[i]],
                   chain = as.factor(chain_id),
                   sub_index = ordinal_suffix(sel_u[i])
        )
      }))
      df_u$parameter <- paste0("u (", df_u$sub_index, " control unit)")
      return(df_u[, c("value", "chain", "parameter")])
    })
    df_u_all <- do.call(rbind, chain_data)
    density_list[["u.sample"]] <- df_u_all
  }
  
  #### 3. f.sample: 선택한 열만 추출 (모든 iteration에서 thin 적용)
  if ("f.sample" %in% names(results[[1]])) {
    chain_data <- lapply(seq_along(results), function(chain_id) {
      mat <- results[[chain_id]][["f.sample"]]
      available_f <- ncol(mat)
      sel_f <- selected_f[selected_f <= available_f]
      thin_idx <- seq(1, nrow(mat), by = thin)
      
      df_f <- do.call(rbind, lapply(seq_along(sel_f), function(i) {
        data.frame(value = mat[thin_idx, sel_f[i]],
                   chain = as.factor(chain_id),
                   sub_index = as.character(sel_f[i])
        )
      }))
      df_f$parameter <- paste0("f (t=", df_f$sub_index, ")")
      return(df_f[, c("value", "chain", "parameter")])
    })
    df_f_all <- do.call(rbind, chain_data)
    density_list[["f.sample"]] <- df_f_all
  }
  
  #### 4. 모든 데이터를 결합하여 facet_wrap으로 그리기
  all_density <- do.call(rbind, density_list)
  
  # --- (A) 부분 문자열 치환: 스칼라 파라미터 레이블 변경
  all_density$parameter <- gsub("sig.sample", "sigma^2", all_density$parameter)
  all_density$parameter <- gsub("tau.sample", "tau", all_density$parameter)
  all_density$parameter <- gsub("l.sample", "kappa", all_density$parameter)
  all_density$parameter <- gsub("w.sample", "eta", all_density$parameter)
  all_density$parameter <- gsub("a.sample", "alpha", all_density$parameter)
  
  # --- (B) 최종 facet 순서 지정 ---
  # 원하는 순서: 스칼라 파라미터(알파, 에타, 시그마^2) → u 샘플 → f 샘플
  final_levels <- c(
    "f (t=1)", "f (t=10)", "f (t=25)", "f (t=50)",
    "tau","eta","alpha", "sigma^2",
    "kappa",
    "u (1st control unit)", "u (10th control unit)", "u (2nd control unit)", "u (7th control unit)"
    
  )
  final_levels <- intersect(final_levels, unique(all_density$parameter))
  all_density$parameter <- factor(all_density$parameter, levels = final_levels)
  
  n_facets <- length(unique(all_density$parameter))
  n_rows <- ceiling(n_facets / 3)
  total_height <- n_rows * 2
  
  p <- ggplot(all_density, aes(x = value, fill = chain)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~ parameter, scales = "free", ncol = 3) +
    labs(title = "",
         x = "Value",
         y = "Density") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      strip.text = element_text(size = 18, face = "bold")
    )
  
  filename <- file.path(output_dir, "density_plots_realdata_0529")
  
  # PNG
  ggsave(paste0(filename, ".png"), plot = p, width = 12, height = total_height, dpi = 150)
  
  # PDF
  ggsave(paste0(filename, ".pdf"), plot = p, width = 12, height = total_height)
  
  # EPS (권장: cairo_ps 사용 → 글꼴/투명도 문제 방지)
  ggsave(paste0(filename, ".eps"), plot = p, width = 12, height = total_height,
         device = cairo_ps, fallback_resolution = 600)
  
  return(p)
}

plot_all_densities(results, thin = 1, 
                   selected_u = c(1, 10, 2, 7), 
                   selected_f = c(1, 10, 25, 50))

#####################################
round(colMeans(beta.sample),2)
round(apply(beta.sample, 2, sd),2)


#####################################################################
####################gamma.sample 시각화##############################
#####################################################################
countries <- c("USA", "UK", "Austria", "Belgium", "Denmark", "France", 
               "Italy", "Netherlands", "Norway", "Switzerland", "Japan", 
               "Greece", "Portugal", "Spain", "Australia", "New Zealand")

barplot(colMeans(results[[2]]$gamma.sample), 
        names.arg = countries,
        main = "Column Means",  # 그래프 제목
        ylab = "",          # y축 라벨
        xlab = "donor pool",
        col = "skyblue")        # 막대 색상 지정



# ggplot2 패키지 설치 및 로드 (설치 안되어 있다면 먼저 설치)
# install.packages("ggplot2")
library(ggplot2)

# 데이터 프레임 생성 (국가 순서를 그대로 유지하기 위해 factor 수준 지정)
df <- data.frame(Country = factor(countries, levels = countries),
                 Mean = colMeans(results[[2]]$gamma.sample))

# ggplot2를 사용한 세로 막대그래프 그리기
ggplot(df, aes(x = Country, y = Mean)) +
  geom_col(fill = "skyblue") +
  labs(title = "", x = "donor pool", y = "gamma") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ggplot2 라이브러리 로드
library(ggplot2)

# 데이터 프레임 생성 (이전에 정의한 countries와 values 벡터를 사용)
df <- data.frame(Country = factor(countries, levels = countries),
                 Mean = values)

# 그래프 생성 (전체 텍스트 크기를 base_size로 16 설정)
plot <- ggplot(df, aes(x = Country, y = Mean)) +
  geom_col(fill = "skyblue") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  labs(title = "", x = "donor pool", y = expression(hat(gamma))) +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# 그래프를 가로로 긴 형식으로 저장 (여기서는 12x6 인치로 저장)
ggsave(filename = "/Users/seul/Desktop/2024/BSC/realdata/260518 revision/plot/realdata_gamma_0529.png", plot = plot, width = 12, height = 4, dpi = 150)
ggsave(filename = "/Users/seul/Desktop/2024/BSC/realdata/260518 revision/plot/realdata_gamma_0529.eps", plot = plot, width = 12, height = 4)
ggsave(filename = "/Users/seul/Desktop/2024/BSC/realdata/260518 revision/plot/realdata_gamma_0529.pdf", plot = plot, width = 12, height = 4)
