#' Geração de ocorrências BARMA
#'
#' Funcao para gerar ocorrências provindas do modelo BARMA com função de ligação arbitrária (aranda, logit, probit, cloglog)
#'
#' @param n Sample size
#' @param freq Number of observations per unit time. E.g. 12 for monthly data
#' @param phi Autoregressive parameter
#' @param theta Moving-average parameter
#' @param beta Covariates
#' @param X Observations of covariates
#' @param alpha Intercept
#' @param prec Precision parameter
#' @param link Link function
#' @param lambda Aranda-Ordaz parameter
#'
#' @export

# n: tamanho da amostra
# freq: frequencia anual de observacoes (12 para observacoes mensais)
simu.barma <- function(n,phi=NA,theta=NA,beta=NA,X=NA,
                            alpha=0.0,prec=NA,freq=12,link="aoz",lambda=NA)
{
  if(is.na(prec)) stop(paste("Precision parameter (prec) must be defined"))

  ar<-NA
  ma<-NA

  if(any(is.na(phi)==F))
  {
    ar <- 1:length(phi)
  }

  if(any(is.na(theta)==F))
  {
    ma <- 1:length(theta)
  }

  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }

  if(linktemp == "aoz"){

    if(is.na(lambda)) stop(paste("Aranda-Ordaz parameter (lambda) must be defined"))
    # source("aux_aoz.r")
    {
      linkfun.aoz <- function(mu,lambda) {
        ret<-log((((1-mu)^(-lambda))-1)/lambda)
        as.vector(ret)
      }

      linkinv.aoz <- function(eta,lambda) {

        eta<- pmax(-.Machine$double.xmax, eta)
        eta<- pmin(.Machine$double.xmax, eta)

        ret<-(1-(1+lambda*exp(eta))^(-1/lambda))

        #print(c(eta,lambda))
        #print(ret)

        ret<- pmin(ret, 1-.Machine$double.eps )
        ret<- pmax(ret, .Machine$double.eps )
        return(as.vector(ret))
      }

      mu.eta.aoz <- function(eta,lambda) {

        eta<- pmax(-.Machine$double.xmax, eta)
        eta<- pmin(.Machine$double.xmax, eta)

        ret<- exp(eta)/((1+lambda*exp(eta))^((1+lambda)/lambda))

        ret<- pmax(ret, .Machine$double.eps)
        as.vector(ret)
      }

      diflink.aoz <- function(y,lambda) {
        ret<- (log(1-y)/(((1-y)^lambda)-1))-(1/lambda)
        as.vector(ret)
      }
    }


    link1 <- structure(list(link = linktemp,
                            linkfun = linkfun.aoz,
                            linkinv = linkinv.aoz
                            # mu.eta = stats$mu.eta.aoz,
                            # diflink = function(t) 1/(stats$mu.eta.aoz(stats$linkfun(t)))
    )
    )

    # linkfun <- link$linkfun
    # linkinv <- link$linkinv
    ###### ARMA model
    if(any(is.na(phi)==F) && any(is.na(theta)==F) && any(is.na(beta)==T))
    {
      #print("ARMA model")
      # print("simu ARMA aoz")

      p <- max(ar)
      q <- max(ma)
      m <- 2*max(p,q)

      ynew <-rep(alpha,(n+m))
      mu <- linkinv.aoz(ynew,lambda)

      error<-rep(0,n+m) # E(error)=0
      eta<- y <- rep(NA,n+m)

      for(i in (m+1):(n+m))
      {
        eta[i]  <- alpha + as.numeric(phi%*%ynew[i-ar]) + as.numeric(theta%*%error[i-ma])
        mu[i]   <- linkinv.aoz(eta[i],lambda)
        y[i]    <- rbeta(1, mu[i]*prec, (1-mu[i])*prec)
        ynew[i] <- linkfun.aoz(y[i],lambda)
        error[i]<- ynew[i]-eta[i]
        #error[i]<- y[i]-mu[i

      }
      return( ts(y[(m+1):(n+m)],frequency=freq) )
    } # ARMA model


    ###### ARMAX model
    if(any(is.na(phi)==F) && any(is.na(theta)==F) && any(is.na(beta)==F))
    {
      # print("simu ARMAX AOZ")
      # seasonal part

      p <- max(ar)
      q <- max(ma)
      m <- 2*max(p,q)
      # m <- max(p,q)


      ynew <-rep(alpha,(n+m))
      mu <- linkinv.aoz(ynew,lambda)

      error<-rep(0,n+m) # E(error)=0
      eta<- y <- rep(NA,n+m)

      X<-as.matrix(X)
      X<-rbind(matrix(0,ncol=dim(X)[2],nrow=m),X)

      for(i in (m+1):(n+m))
      {
        eta[i]  <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) + (theta%*%error[i-ma])
        error[i] <- ynew[i]-eta[i]  #ma - theta
        mu[i]   <- linkinv.aoz(eta[i],lambda)
        y[i]    <- rbeta(1, mu[i]*prec, (1-mu[i])*prec)
        ynew[i] <- linkfun.aoz(y[i],lambda)
        error[i]<- ynew[i]-eta[i]
      }
      # print("ARMAX")
      return(ts(y[(m+1):(n+m)],frequency=freq))
    } # ARMAX model


    ###### AR model
    if(any(is.na(phi)==F) && any(is.na(theta)==T) && any(is.na(beta)==T))
    {
      #print("AR model")

      p <- max(ar)
      m <- 2*p

      ynew <-rep(alpha,(n+m))
      mu <- linkinv.aoz(ynew,lambda)

      eta<- y <- rep(NA,n+m)

      for(i in (m+1):(n+m))
      {
        eta[i]  <- alpha + (phi%*%ynew[i-ar])
        mu[i]   <- linkinv.aoz(eta[i],lambda)
        y[i]    <- rbeta(1, mu[i]*prec, (1-mu[i])*prec)
        ynew[i] <- linkfun.aoz(y[i],lambda)
      }
      print("AR")
      return( ts(y[(m+1):(n+m)],frequency=freq) )
    } # AR model


    ###### MA model
    if(any(is.na(phi)==T) && any(is.na(theta)==F) && any(is.na(beta)==T))
    {
      #print("MA model")

      q <- max(ma)
      m <- 2*q

      ynew <-rep(alpha,(n+m))
      mu <- linkinv.aoz(ynew,lambda)

      eta <- y <- error <- rep(0,n+m) # E(error)=0

      #print(ma)

      for(i in (m+1):(n+m))
      {
        eta[i]  <- alpha + (theta%*%error[i-ma])
        mu[i]   <- linkinv.aoz(eta[i],lambda)
        y[i]    <- rbeta(1, mu[i]*prec, (1-mu[i])*prec)
        ynew[i] <- linkfun.aoz(y[i],lambda)
        error[i]<- ynew[i]-eta[i]
        #error[i]<- y[i]-mu[i]
      }
      print("MA")
      return( ts(y[(m+1):(n+m)],frequency=freq) )
    } # fim MA model

  }else{

    # linktemp <- substitute(link)
    # if (!is.character(linktemp))
    # {
    #   linktemp <- deparse(linktemp)
    #   if (linktemp == "link")
    #     linktemp <- eval(link)
    # }
    if (any(linktemp == c("logit", "probit", "cloglog")))
    {
      stats <- make.link(linktemp)
    }else{
      stop(paste(linktemp, "link not available, available links are \"logit\", ",
                 "\"probit\" and \"cloglog\""))
    }

    link <- structure(list(link = linktemp,
                           linkfun = stats$linkfun,
                           linkinv = stats$linkinv
    )
    )

    linkfun <- link$linkfun
    linkinv <- link$linkinv
  }

  ###### ARMA model
  if(any(is.na(phi)==F) && any(is.na(theta)==F))
  {
    print("ARMA logit")
    # seasonal part

    p <- max(ar)
    q <- max(ma)
    m <- 2*max(p,q)

    ynew <-rep(alpha,(n+m))
    mu <- linkinv(ynew)

    error<-rep(0,n+m) # E(error)=0
    eta<- y <- rep(NA,n+m)

    for(i in (m+1):(n+m))
    {
      eta[i]  <- alpha + as.numeric(phi%*%ynew[i-ar]) + as.numeric(theta%*%error[i-ma])
      mu[i]   <- linkinv(eta[i])
      y[i]    <- rbeta(1, mu[i]*prec, (1-mu[i])*prec)
      ynew[i] <- linkfun(y[i])
      error[i]<- ynew[i]-eta[i]
      #error[i]<- y[i]-mu[i]
    }

    return( ts(y[(m+1):(n+m)],frequency=freq) )
  } # ARMA model


  ###### AR model
  if(any(is.na(phi)==F) && any(is.na(theta)==T))
  {
    print("AR model")

    p <- max(ar)
    m <- 2*p

    ynew <-rep(alpha,(n+m))
    mu <- linkinv(ynew)

    eta<- y <- rep(NA,n+m)

    for(i in (m+1):(n+m))
    {
      eta[i]  <- alpha + (phi%*%ynew[i-ar])
      mu[i]   <- linkinv(eta[i])
      y[i]    <- rbeta(1, mu[i]*prec, (1-mu[i])*prec)
      ynew[i] <- linkfun(y[i])
    }

    return( ts(y[(m+1):(n+m)],frequency=freq) )
  } # AR model


  ###### MA model
  if(any(is.na(phi)==T) && any(is.na(theta)==F))
  {
    print("MA model")

    q <- max(ma)
    m <- 2*q

    ynew <-rep(alpha,(n+m))
    mu <- linkinv(ynew)

    eta <- y <- error <- rep(0,n+m) # E(error)=0

    #print(ma)

    for(i in (m+1):(n+m))
    {
      eta[i]  <- alpha + (theta%*%error[i-ma])
      mu[i]   <- linkinv(eta[i])
      y[i]    <- rbeta(1, mu[i]*prec, (1-mu[i])*prec)
      ynew[i] <- linkfun(y[i])
      error[i]<- ynew[i]-eta[i]
      #error[i]<- y[i]-mu[i]
    }

    return( ts(y[(m+1):(n+m)],frequency=freq) )
  } # fim MA model

}
#
# set.seed(369)
# simu.barma.full(60,phi=c(0.45,0.3),theta=c(0.5),prec=100)
# set.seed(369)
# simu.barma.full(60,phi=c(0.45,0.3),theta=c(0.5),prec=100,link = "logit") #alpha=0.0
# set.seed(369)
# simu.barma.full(60,phi=c(0.45,0.3),theta=c(0.5),prec=100,link = "aoz",lambda=1.5) #alpha


