# Implementado por Fabio M Bayer (bayer@ufsm.br) em 02/10/2014
# Alteracoes em 17/06/2015
#
# Algumas informacoes:
# diag = 0 : nao plota graficos
# diag = 1 : plota graficos na tela
# diga = 2 : gera graficos em pdf e plota graficos na tela
#
# h : quantidade de passos a frente para fazer previsao
#
# O objeto de entrada da funcao deve ser serie temporal (ts)
#
# Exemplos de uso:
#
# 1) BARMA(2,3) com funcao de ligacao logit e previsao 6 passos a frente
# fit <- barma(y,ar=c(1,2),ma=c(1,2))
# 
# 2) SBARMA(p,q)(P,Q) com funcao de ligacao cloglog e previsao 12 passos a frente
# fit <- barma(y,ar=c(1,2,5),ma=c(1,3),AR=c(1,2),MA=1,link="cloglog",h=12)
# 
# Obs: Perceba que pode ser os lags que voce desejar.
# 
# 3) imprimindo graficos em arquivos pdf
# fit <- barma(y,ar=c(1,2),ma=c(1,2),diag=2)


barma<- function (y, ar=NA, ma=NA, AR=NA, MA=NA, link = "logit",diag=1,h=6,X=NA,X_hat=NA)
{  
  source("barma.fit.r")
  
  if (min(y) <= 0 || max(y) >= 1)
    stop("OUT OF RANGE (0,1)!")
  
  if(is.ts(y)==T)
  {
    freq<-frequency(y)
  }else stop("data can be a time-series object")
  
  
  if(any(is.na(ar))==F)
  {
    if(any(is.na(AR)))
    {
      names_phi<-c(paste("phi",ar,sep=""))
    }else{
      names_phi<-c(paste("phi",ar,sep=""),paste("PHI",AR,sep=""))
      ar<- c(ar,freq*AR)
    }
  }else{
    if(any(is.na(AR)))
    {
      names_phi <- NA
    }else{
      names_phi<-c(paste("PHI",AR,sep=""))
      ar<- c(freq*AR)
    }
  }
  
  if(any(is.na(ma))==F)
  {
    if(any(is.na(MA)))
    {
      names_theta<-c(paste("theta",ma,sep=""))
    }else{
      names_theta<-c(paste("theta",ma,sep=""),paste("THETA",MA,sep=""))
      ma<- c(ma,freq*MA)
    }
  }else{
    if(any(is.na(MA)))
    {
      names_theta <- NA
    }else{
      names_theta<-c(paste("THETA",MA,sep=""))
      ma<- c(freq*MA)
    }
  }
  
  if(any(is.na(X))==F)
  {
    names_beta<-c(paste("beta",1:ncol( as.matrix(X) ),sep=""))
  }
  
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  
  m <- max(p,q,na.rm=T)
  
  p1 <- length(ar)
  q1 <- length(ma)
  
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("logit", "probit", "cloglog")))
    stats <- make.link(linktemp)
  else stop(paste(linktemp, "link not available, available links are \"logit\", ",
                  "\"probit\" and \"cloglog\""))
  
  link1 <- structure(list(link = linktemp, 
                          linkfun = stats$linkfun,
                          linkinv = stats$linkinv, 
                          mu.eta = stats$mu.eta, 
                          diflink = function(t) 1/(stats$mu.eta(stats$linkfun(t)))
  )
  )
  
    
  fit1 <- barma.fit(y, ar, ma, link1, names_phi, names_theta, names_beta, diag, h, X, X_hat) # model estimation
  
  return(fit1)
}

