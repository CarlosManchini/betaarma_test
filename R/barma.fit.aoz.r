# SCRIPT AQUELE modified

barma.fit.aoz<- function (y, ar, ma, link, names_phi,names_theta,names_beta,diag,h1,X,X_hat,resid,lambda=NA)
{
  maxit1<-1000
  z<-c()
  cont_conv<-0

  source("aux_aoz.r")
  
  linkfun <- linkfun.aoz  #changed
  linkinv <- linkinv.aoz
  mu.eta <-  mu.eta.aoz
  diflink <- diflink.aoz
  mu.eta <- mu.eta.aoz
  
  ynew = linkfun(y,lambda)
  ystar = log(y/(1-y))
  
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  m <- max(p,q,na.rm=T)
  p1 <- length(ar)
  q1 <- length(ma)
  
  y_prev <- c(rep(NA,(n+h1)))
  eta<-rep(NA,n)
  
  # inicializacao dos parametros alpha e phi (beta)
  if(any(is.na(ar)==F)) # se TEM componente ar
  {
    P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
    
    for(i in 1:(n-m))
    {
      P[i,] <- ynew[i+m-ar]
    }
    
    Z <- cbind(rep(1,(n-m)),P)
  }else{
    Z <- as.matrix(rep(1,(n-m)))
  }
  
  if(any(is.na(X)==T)) # nao tem regressores
  {
    x <- as.matrix(Z)
    Y <- y[(m+1):n]
    Ynew = linkfun(Y,lambda)
    Ystar = log(Y/(1-Y))
    ajuste = lm.fit(x, Ynew)
    mqo = c(ajuste$coef)
    k = length(mqo)
    n1 = length(Y)
    meangy = fitted(ajuste)
    meany = linkinv(meangy,lambda)
    er<-Y-meany
    vary<-var(er)
    prec = mean((meany*(1-meany))/vary-1)
  }else{ # com regressores
    X_hat<-as.matrix(X_hat)
    X<-as.matrix(X)
    x <- cbind(as.matrix(Z),X[(m+1):n,])
    Y <- y[(m+1):n]
    Ynew = linkfun(Y,lambda)
    Ystar = log(Y/(1-Y))
    ajuste = lm.fit(x, Ynew)
    mqo = c(ajuste$coef)
    k = length(mqo)
    n1 = length(Y)
    meangy = fitted(ajuste)
    meany = linkinv(meangy,lambda)
    er<-Y-meany
    vary<-var(er)
    prec = mean((meany*(1-meany))/vary-1)
    # print(prec)
  }
  
  ############
  
  ######### BARMA X model
  if(any(is.na(ar)==F) && any(is.na(ma)==F) && any(is.na(X)==F))
  { 
     print("BARMA~X model aoz") #dps de conferir, comentar
    beta1<- mqo[(p1+2):length(mqo)]
    reg <- c(mqo[1:(p1+1)], rep(0,q1), prec, lambda ,beta1) # 2*prec  #initializing the parameter values
    
    loglik <- function(z) 
    {
      alpha <- z[1]
      phi = z[2:(p1+1)] 
      theta = z[(p1+2):(p1+q1+1)]
      prec <- z[p1+q1+2] # precision parameter
      lambda <- z[p1+q1+3]  ##added
      beta <- z[(p1+q1+4):length(z)]
      
      error<-rep(0,n) # E(error)=0 
      eta<-rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) + (theta%*%error[i-ma])
        error[i] <- ynew[i]-eta[i]  #ma - theta
        #error[i] <- y[i]-linkinv(eta[i])
        #mui<-linkinv(eta[i])
        #error[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
      }
      mu <- linkinv(eta[(m+1):n],lambda)
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings(dbeta(y1, mu * prec, (1 - mu) * prec, log = TRUE))
      sum(ll)
    } 
    
    {
      escore <- function(z)
      {
        alpha <- z[1]
        phi = z[2:(p1+1)]
        theta = z[(p1+2):(p1+q1+1)]
        prec <- z[p1+q1+2] # precision parameter
        lambda <- z[p1+q1+3]  ##added
        beta <- z[(p1+q1+4):length(z)]

        error<- rep(0,n) # E(error)=0
        eta<- rep(NA,n)

        for(i in (m+1):n)
        {
          eta[i]<- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) + (theta%*%error[i-ma])
          error[i]<- ynew[i]-eta[i] #MA - theta
          #error[i] <- y[i]-linkinv(eta[i])
          #mui<-linkinv(eta[i])
          #error[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
        }

        mu <- linkinv(eta[(m+1):n],lambda)
        y1<-y[(m+1):n]
        ystar <- log(y1/(1-y1)) #
        mustar <- digamma(mu * prec) - digamma((1 - mu) * prec) #

        mT <- diag(mu.eta(eta[(m+1):n],lambda)) # corrigir (mu) por (eta[(m+1):n])

        R <- matrix(rep(NA,(n-m)*q1),ncol=q1) # MA - theta
        for(i in 1:(n-m))
        {
          R[i,] <- error[i+m-ma] #P[i,]:AR feito no topo
        }

        P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
        for(i in 1:(n-m))
        {
          P[i,] <- ynew[i+m-ar] - X[i+m-ar,]%*%as.matrix(beta)
        }

        k1<- length(beta)
        M <- matrix(rep(NA,(n-m)*k1),ncol=k1)
        for(i in 1:(n-m))
        {
          for(j in 1:k1) #comentado no karma
            M[i,j] <- X[i+m,j]-sum(phi*X[i+m-ar,j]) #comentado no karma

          # M[i,] <- X[i+m,]-(phi%*%X[i+m-ar,]) #usado no karma              # VER

        }

        ###FB  recorrencias
        deta.dalpha <- rep(0,n)
        deta.dphi <- matrix(0, ncol=p1,nrow=n)
        deta.dtheta<- matrix(0, ncol=q1,nrow=n)
        deta.dlambda<- rep(0,n) #matrix(0, ncol=1, nrow=n) #me
        dmu.dlambda<- rep(0,n)
        deta.dbeta<- matrix(0, ncol=k1,nrow=n)

        # print(dim(P))
        # print(dim(deta.dphi))

        for(i in (m+1):n)
        {
          deta.dalpha[i] <- 1 - theta%*%deta.dalpha[i-ma]
          deta.dphi[i,] <- P[(i-m),] - theta%*%deta.dphi[i-ma,]
          deta.dtheta[i,] <- R[(i-m),] - theta%*%deta.dtheta[i-ma,]

          dmu.dlambda[i] <- (((1+lambda*exp(eta[i]))^(-1/lambda))/lambda)*( (1/(exp(-eta[i])+lambda))*(1+deta.dlambda[i])-(log(1+lambda*exp(eta[i]))/lambda) )

          deta.dbeta[i,]<- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
        }

        a <- deta.dalpha[(m+1):n]
        rP <- deta.dphi[(m+1):n,]
        rR <- deta.dtheta[(m+1):n,]
        #rL <- deta.dlambda[(m+1):n]
        rL <- dmu.dlambda[(m+1):n]
        rM <- deta.dbeta[(m+1):n,]

        Ualpha <- prec * a %*% mT %*% (ystar-mustar)
        Uphi <-   prec * t(rP) %*% mT %*% (ystar-mustar)
        Utheta <- prec * t(rR) %*% mT %*% (ystar-mustar)
        Uprec <-  sum(mu * (ystar-mustar) + log(1-y1)
                      - digamma((1 - mu) * prec) + digamma(prec) )
        # Ulambda <- prec * rL %*% mT %*% (ystar-mustar)
        Ulambda <- prec *sum((ystar-mustar) * rL)
        Ubeta <-   prec * t(rM) %*% mT %*% (ystar-mustar)

        rval <- c(Ualpha,Uphi,Utheta,Uprec,Ulambda,Ubeta)
      }
    } #escore..
    
    names_par <- c("alpha",names_phi,names_theta,"precision","lambda",names_beta)
    
    # print(reg) # intercepto, 2 autoregressivos (phi), medias moveis (theta), precisao e lambda
    
    opt <- optim(reg, loglik, escore,
                 method = "BFGS",
                 control = list(fnscale = -1, maxit = 1000000, reltol = 1e-12))
    
    # opt <- optim(reg, loglik, escore,
    #              method = "L-BFGS-B",
    #              lower = c(rep(-Inf,p1+q1+1),0.01,0.01,rep(-Inf,dim(X)[2])),
    #              upper = rep(Inf,p1+q1+3+dim(X)[2]),
    #              control = list(fnscale = -1, maxit = 100000))#, reltol = 1e-12))
    
    if (opt$conv != 0)
    {
      warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
      opt <- optim(reg, loglik, #escore,
                   method = "BFGS",
                   control = list(fnscale = -1, maxit = 10000000, reltol = 1e-12))
      
      # opt <- optim(reg, loglik, #escore,
      #              method = "L-BFGS-B",
      #              lower = c(rep(-Inf,p1+q1+1),0.1,0.1,rep(-Inf,dim(X)[2])),
      #              upper = rep(Inf,p1+q1+3+dim(X)[2]),
      #              control = list(fnscale = -1, maxit = 100000))#, reltol = 1e-12))
      
      cont_conv <- cont_conv+1
    }
    
    # library(rootSolve)
    # escores<-rbind(escore(opt$par),gradient(loglik, opt$par))
    # colnames(escores)<-names_par
    # print(escores) # deve ser prox de zero
    
    # print("AQUI6")
    # print(opt$par)
    
    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+q1+3+ncol(X))]  #[1:(p1+q1+2)]
    names(coef)<-names_par
    z$coeff <- coef
    z$contador <- cont_conv
    
    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    theta <- coef[(p1+2):(p1+q1+1)]
    prec <- coef[p1+q1+2] # precision parameter
    lambda <- coef[p1+q1+3]
    beta <- coef[(p1+q1+4):length(coef)]
    
    
    z$alpha <- alpha
    z$phi <- phi
    z$theta <- theta
    z$prec <- prec
    z$lambda <- lambda
    z$beta <- beta
    
   
    errorhat<-rep(0,n) # E(error)=0
    etahat<-rep(NA,n)

    for(i in (m+1):n)
    {
      etahat[i]<-alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) + (theta%*%errorhat[i-ma])
      errorhat[i]<- ynew[i]-etahat[i] # predictor scale # ma -theta

      #errorhat[i] <- y[i]-linkinv(etahat[i]) # original scale
      #mui<-linkinv(etahat[i])
      #errorhat[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
    }
    muhat <- linkinv(etahat[(m+1):n],lambda)
    y1<-y[(m+1):n]

    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    z$mustarhat <- digamma(muhat * prec) - digamma((1 - muhat) * prec)

    # inicio fisher
    {
    R <- matrix(rep(NA,(n-m)*q1),ncol=q1) # MA - theta
    for(i in 1:(n-m))
    {
      R[i,] <- errorhat[i+m-ma]
    }

    P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
    for(i in 1:(n-m))
    {
      P[i,] <- ynew[i+m-ar] - X[i+m-ar,]%*%as.matrix(beta)
    }

    k1<-length(beta)
    M <- matrix(rep(NA,(n-m)*k1),ncol=k1)
    for(i in 1:(n-m))
    {
      for(j in 1:k1) #aqui NAO eh comentado no karma
        M[i,j] <- X[i+m,j]-sum(phi*X[i+m-ar,j])
    }

    vI <- as.vector(rep(1,n-m)) #wtff are this here?



    ###FB  recorrencias  #FISHER
    deta.dalpha <- rep(0,n)
    deta.dphi <- matrix(0, ncol=p1,nrow=n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)
    deta.dlambda<- rep(0,n) #matrix(0, ncol=1, nrow=n) #me
    dmu.dlambda<- rep(0,n)
    deta.dbeta<- matrix(0, ncol=k1,nrow=n)

    # print(dim(P))
    # print(dim(deta.dphi))

    for(i in (m+1):n) #FISHER
    {
      deta.dalpha[i] <- 1 - theta%*%deta.dalpha[i-ma]
      deta.dphi[i,] <- P[(i-m),] - theta%*%deta.dphi[i-ma,]
      deta.dtheta[i,] <- R[(i-m),] - theta%*%deta.dtheta[i-ma,]

      dmu.dlambda[i] <- (((1+lambda*exp(etahat[i]))^(-1/lambda))/lambda)*( (1/(exp(-etahat[i])+lambda))*(1+deta.dlambda[i])-(log(1+lambda*exp(etahat[i]))/lambda) )

      deta.dbeta[i,]<- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
    }

    a <- deta.dalpha[(m+1):n]
    rP <- deta.dphi[(m+1):n,]
    rR <- deta.dtheta[(m+1):n,]
    #rL <- deta.dlambda[(m+1):n]
    rL <- dmu.dlambda[(m+1):n]
    rM <- deta.dbeta[(m+1):n,]

    psi1 = trigamma(muhat * prec)
    psi2 = trigamma((1 - muhat) * prec)
    mT <- diag(mu.eta( etahat[(m+1):n] ,lambda)) #muhat por etahat[(m+1):n]
   # W = diag(c(prec * (psi1 + psi2))) %*% mT^2   #v_t * T^2
    VT2 = diag(c(prec^2 * (psi1 + psi2))) %*% mT^2   #v_t * T^2
    VT = diag(c(prec^2 * (psi1 + psi2))) %*% mT   #v_t * T
    V = diag(c(prec^2 * (psi1 + psi2)))   #v_t
    vc = prec * (psi1 * muhat - psi2 * (1 - muhat))  #c_t
    D = diag(c(psi1 * (muhat^2) + psi2 * (1 - muhat)^2 - trigamma(prec))) # S


    Kaa <-  t(a) %*% VT2 %*% a
    Kpa <-  t(rP) %*% VT2 %*% a #  t(a) %*% VT2 %*% rP
    Kap <- t(Kpa)
    Kta <-  t(rR) %*% VT2 %*% a # t(a) %*% VT2 %*% rR
    Kat <- t(Kta)
    Kaprec <- t(a) %*% mT %*% vc
    Kpreca <- Kaprec
    Kpp <-  t(rP) %*% VT2 %*% rP
    Kpt <-  t(rP) %*% VT2 %*% rR
    Ktp <- t(Kpt)
    Kpprec <- t(rP) %*% mT %*% vc
    Kprecp <- t(Kpprec)
    Ktt <- t(rR) %*% VT2 %*% rR
    Ktprec <- t(rR) %*% mT %*% vc
    Kprect <- t(Ktprec)
    Kprecprec <- sum(diag(D))


    Kba <- t(rM)%*% VT2 %*% a # t(a) %*% VT2 %*% rM
    Kbb <-  t(rM) %*% VT2 %*% rM
    Kbprec <- t(rM) %*% mT %*% vc
    Kbp <- t(rM) %*% VT2 %*% rP
    Kbt <- t(rM) %*% VT2 %*% rR

    Kab <- t(Kba)
    Kprecb <- t(Kbprec)
    Kpb <- t(Kbp)
    Ktb <- t(Kbt)

    Kal <- t(a) %*% VT %*% rL
    Kla <- t(Kal)
    Kbl <- t(rM) %*% VT %*% rL
    Klb <- t(Kbl)
    Kpl <- t(rP) %*% VT %*% rL
    Klp <- t(Kpl)
    Ktl <- t(rR) %*% VT %*% rL #rP nu fodi
    Klt <- t(Ktl)
    Kprecl <- vc %*% rL
    Klprec <- t(Kprecl)
    Kll <- t(rL) %*% V %*% rL

    K <- rbind(
      cbind(Kaa,Kap,Kat,Kaprec,Kal,Kab),
      cbind(Kpa,Kpp,Kpt,Kpprec,Kpl,Kpb),
      cbind(Kta,Ktp,Ktt,Ktprec,Ktl,Ktb),
      cbind(Kpreca,Kprecp,Kprect,Kprecprec,Kprecl,Kprecb),
      cbind(Kla,Klp,Klt,Klprec,Kll,Klb),
      cbind(Kba,Kbp,Kbt,Kbprec,Kbl,Kbb)
          )

    z$fisher <- K

    }
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted

    X_prev<- rbind(X,X_hat)

    #     (phi%*%(ynew[i-ar]-X[i-ar,]%*%beta ))
    #     (phi%*%(ynew_prev[n+i-ar]-X_prev[i-ar,]%*%beta ))

    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) + (phi%*%(ynew_prev[n+i-ar]-X_prev[n+i-ar,]%*%as.matrix(beta) )) + (theta%*%errorhat[n+i-ma])
      #errorhat[i]<-ynew[i]-etahat[i] # predictor scale
      y_prev[n+i] <- linkinv(ynew_prev[n+i],lambda)
      errorhat[n+i] <- 0 # original scale
    }
  } 
  
  ######### BARMA model
  if(any(is.na(ar)==F) && any(is.na(ma)==F) && any(is.na(X)==T))
  { 
    # print("BARMA model aoz",quote=F)
    
    reg <- c(mqo, rep(0,q1), prec, lambda) # initializing the parameter values
    
    loglik <- function(z) 
    {
      alpha <- z[1]
      phi = z[2:(p1+1)] 
      theta = z[(p1+2):(p1+q1+1)]
      prec <- z[p1+q1+2] # precision parameter
      lambda <- z[p1+q1+3]
      
      error<-rep(0,n) # E(error)=0 
      eta<-rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i] <- alpha + (phi%*%ynew[i-ar]) + (theta%*%error[i-ma])
        error[i] <- ynew[i]-eta[i] 
        #error[i] <- y[i]-linkinv(eta[i])
        #mui<-linkinv(eta[i])
        #error[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
      }
      mu <- linkinv(eta[(m+1):n], lambda)
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings(dbeta(y1, mu * prec, (1 - mu) * prec, log = TRUE))
      sum(ll)
    } 
    
    escore <- function(z) 
    {
      alpha <- z[1]
      phi <- z[2:(p1+1)]
      theta <- z[(p1+2):(p1+q1+1)]
      prec <- z[p1+q1+2] # precision parameter
      lambda <- z[p1+q1+3]
      
      
      error<- rep(0,n) # E(error)=0 
      eta<- rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i]<- alpha + (phi%*%ynew[i-ar]) + (theta%*%error[i-ma])
        error[i]<- ynew[i]-eta[i] 
        #error[i] <- y[i]-linkinv(eta[i])
        #mui<-linkinv(eta[i])
        #error[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
      }
      
      mu <- linkinv(eta[(m+1):n], lambda)
      y1<-y[(m+1):n]
      ystar <- log(y1/(1-y1))
      mustar <- digamma(mu * prec) - digamma((1 - mu) * prec)
      
      mT <- diag(mu.eta(eta[(m+1):n], lambda))
      
      R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
      for(i in 1:(n-m))
      {
        R[i,] <- error[i+m-ma]
      }
      
      ###FB  recorrencias
      deta.dalpha <- rep(0,n)
      deta.dphi <- matrix(0, ncol=p1,nrow=n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)
      dmu.dlambda <- rep(0,n)
      deta.dlambda <- rep(0,n)
      # print(dim(P))
      # print(dim(deta.dphi))
      
      for(i in (m+1):n)
      {
        #print(P[(i-m),(i-m)])
        deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
        deta.dphi[i,]<- P[(i-m),] - theta%*%deta.dphi[i-ma,]
        deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
        dmu.dlambda[i] <- (((1+lambda*exp(eta[i]))^(-1/lambda))/lambda)*( (1/(exp(-eta[i])+lambda))*(1+deta.dlambda[i])-(log(1+lambda*exp(eta[i]))/lambda) )
        
      }
      
      a <- deta.dalpha[(m+1):n]
      rP <- deta.dphi[(m+1):n,]
      rR <- deta.dtheta[(m+1):n,]
      rL <- dmu.dlambda[(m+1):n]
      
      
      Ualpha <- prec * a %*% mT %*% (ystar-mustar)
      Uphi <-   prec * t(rP) %*% mT %*% (ystar-mustar)
      Utheta <- prec * t(rR) %*% mT %*% (ystar-mustar)
      Uprec <-  sum(mu * (ystar-mustar) + log(1-y1)
                    - digamma((1 - mu) * prec) + digamma(prec) )
      #Ulambda <- prec * rL %*% mT %*% (ystar-mustar)
      Ulambda <- prec *sum((ystar-mustar) * rL)
      
      rval <- c(Ualpha,Uphi,Utheta,Uprec,Ulambda)
    }
    names_par <- c("alpha", names_phi, names_theta, "precision", "lambda")
    
    opt <- optim(reg, loglik, escore,
                 method = "BFGS",
                 control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
    

    
    if (opt$conv != 0)
    {
      warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
      opt <- optim(reg, loglik, #escore,
                   method = "BFGS",
                   control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))

      
      cont_conv <- cont_conv+1
    }
    
    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+q1+3)]
    names(coef)<-names_par
    z$coeff <- coef
    z$contador <- cont_conv
    
    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    theta <- coef[(p1+2):(p1+q1+1)]
    prec <- coef[p1+q1+2] # precision parameter
    lambda <- coef[p1+q1+3]
    
    
    z$alpha <- alpha
    z$phi <- phi
    z$theta <- theta
    z$prec <- prec
    z$lambda <- lambda
    
    
    errorhat<-rep(0,n) # E(error)=0
    etahat<-rep(NA,n)
    
    for(i in (m+1):n)
    {
      etahat[i]<-alpha + (phi%*%ynew[i-ar]) + (theta%*%errorhat[i-ma])
      errorhat[i]<- ynew[i]-etahat[i] # predictor scale
      #errorhat[i] <- y[i]-linkinv(etahat[i]) # original scale
      #mui<-linkinv(etahat[i])
      #errorhat[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
    }
    muhat <- linkinv(etahat[(m+1):n],lambda) #vou ficar por aqui mesmo hu3hu3
    y1<-y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    z$mustarhat <- digamma(muhat * prec) - digamma((1 - muhat) * prec)
    {
      R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
      for(i in 1:(n-m))
      {
        R[i,] <- errorhat[i+m-ma]
      }
      
      vI <- as.vector(rep(1,n-m))
      # 
      # 
      ###FB  recorrencias  #FISHER
      deta.dalpha <- rep(0,n)
      deta.dphi <- matrix(0, ncol=p1,nrow=n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)
      deta.dlambda<- rep(0,n) #matrix(0, ncol=1, nrow=n) #me
      dmu.dlambda<- rep(0,n)
      # 
      # # print(dim(P))
      # # print(dim(deta.dphi))
      # 
      for(i in (m+1):n) #FISHER
      {
        deta.dalpha[i] <- 1 - theta%*%deta.dalpha[i-ma]
        deta.dphi[i,] <- P[(i-m),] - theta%*%deta.dphi[i-ma,]
        deta.dtheta[i,] <- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
        
        dmu.dlambda[i] <- (((1+lambda*exp(etahat[i]))^(-1/lambda))/lambda)*( (1/(exp(-etahat[i])+lambda))*(1+deta.dlambda[i])-(log(1+lambda*exp(etahat[i]))/lambda) )
      }
      a <- deta.dalpha[(m+1):n]
      rP <- deta.dphi[(m+1):n,]
      rR <- deta.dtheta[(m+1):n,]
      #rL <- deta.dlambda[(m+1):n]
      rL <- dmu.dlambda[(m+1):n]
      
      psi1 = trigamma(muhat * prec)
      psi2 = trigamma((1 - muhat) * prec)
      mT <- diag(mu.eta( etahat[(m+1):n] ,lambda)) #muhat por etahat[(m+1):n]
      # W = diag(c(prec * (psi1 + psi2))) %*% mT^2   #v_t * T^2
      VT2 = diag(c(prec^2 * (psi1 + psi2))) %*% mT^2   #v_t * T^2
      VT = diag(c(prec^2 * (psi1 + psi2))) %*% mT   #v_t * T
      V = diag(c(prec^2 * (psi1 + psi2)))   #v_t
      vc = prec * (psi1 * muhat - psi2 * (1 - muhat))  #c_t
      D = diag(c(psi1 * (muhat^2) + psi2 * (1 - muhat)^2 - trigamma(prec))) # S
      
      
      Kaa <-  t(a) %*% VT2 %*% a
      Kpa <-  t(rP) %*% VT2 %*% a #  t(a) %*% VT2 %*% rP
      Kap <- t(Kpa)
      Kta <-  t(rR) %*% VT2 %*% a # t(a) %*% VT2 %*% rR
      Kat <- t(Kta)
      Kaprec <- t(a) %*% mT %*% vc
      Kpreca <- Kaprec
      Kpp <-  t(rP) %*% VT2 %*% rP
      Kpt <-  t(rP) %*% VT2 %*% rR
      Ktp <- t(Kpt)
      Kpprec <- t(rP) %*% mT %*% vc
      Kprecp <- t(Kpprec)
      Ktt <- t(rR) %*% VT2 %*% rR
      Ktprec <- t(rR) %*% mT %*% vc
      Kprect <- t(Ktprec)
      Kprecprec <- sum(diag(D))
      
      
      Kal <- t(a) %*% VT %*% rL
      Kla <- t(Kal)
      Kpl <- t(rP) %*% VT %*% rL
      Klp <- t(Kpl)
      Ktl <- t(rR) %*% VT %*% rL #rP nu fodi
      Klt <- t(Ktl)
      Kprecl <- vc %*% rL
      Klprec <- t(Kprecl)
      Kll <- t(rL) %*% V %*% rL
      
      K <- rbind(
        cbind(Kaa,Kap,Kat,Kaprec,Kal),
        cbind(Kpa,Kpp,Kpt,Kpprec,Kpl),
        cbind(Kta,Ktp,Ktt,Ktprec,Ktl),
        cbind(Kpreca,Kprecp,Kprect,Kprecprec,Kprecl),
        cbind(Kla,Klp,Klt,Klprec,Kll)
      )
      
      z$fisher <- K
      
      #### Forecasting
      ynew_prev <- c(ynew,rep(NA,h1))
      y_prev[1:n] <- z$fitted
      
      for(i in 1:h1)
      {
        ynew_prev[n+i] <- alpha + (phi%*%ynew_prev[n+i-ar]) + (theta%*%errorhat[n+i-ma])
        #errorhat[i]<-ynew[i]-etahat[i] # predictor scale
        y_prev[n+i] <- linkinv(ynew_prev[n+i],lambda)
        errorhat[n+i] <- 0 # original scale
      }
    } #fisher
    
  }
  
  ##### only BAR model
  if(any(is.na(ar)==F) && any(is.na(ma)==T) && any(is.na(X)==T))
  {
    print("BAR model",quote=F)
    q1<-0
    reg <- c(mqo, prec,lambda) # initializing the parameter values
    
    loglik <- function(z) 
    {
      alpha <- z[1]
      phi = z[2:(p1+1)] 
      prec <- z[p1+2] # precision parameter
      lambda<-z[p1+3]
      
      eta<-rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i]<-alpha + (phi%*%ynew[i-ar]) 
      }
      mu <- linkinv(eta[(m+1):n],lambda)
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings(dbeta(y1, mu * prec, (1 - mu) * prec, log = TRUE))
      sum(ll)
    } 
    
    escore <- function(z) 
    {
      alpha <- z[1]
      phi <- z[2:(p1+1)]
      prec <- z[p1+2] # precision parameter
      lambda<-z[p1+3]
      
      eta<- rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i]<- alpha + (phi%*%ynew[i-ar]) 
      }
      
      mu <- linkinv(eta[(m+1):n],lambda)
      y1<-y[(m+1):n]
      ystar <- log(y1/(1-y1))
      mustar <- digamma(mu * prec) - digamma((1 - mu) * prec)
      
      mT <- diag(mu.eta(eta[(m+1):n],lambda)) #changed
      rL <- (((1+lambda*exp(eta[(m+1):n]))^(-1/lambda))/lambda)*( (1/(exp(-eta[(m+1):n])+lambda))*(1+0)-(log(1+lambda*exp(eta[(m+1):n]))/lambda) ) 

      Ualpha <- prec * sum((mu.eta(eta[(m+1):n],lambda)) * (ystar-mustar))  ##mu por eta[(m+1):n] sera?
      Uphi <-   prec * t(P) %*% mT %*% (ystar-mustar)
      Uprec <-  sum(mu * (ystar-mustar) + log(1-y1) 
                    - digamma((1 - mu) * prec) + digamma(prec) )
      Ulambda <- prec * sum((ystar-mustar) * rL)
      
      rval <- c(Ualpha,Uphi,Uprec,Ulambda)
    }
    names_par <- c("alpha",names_phi,"precision","lambda")
    
    opt <- optim(reg, loglik, escore, method = "BFGS", 
                 control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
    
    if (opt$conv != 0) ###
    {
      warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
      opt <- optim(reg, loglik, #escore,
                   method = "BFGS", 
                   control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
      if (opt$conv != 0)
      {
        warning("FUNCTION DID NOT CONVERGE NEITHER WITH NUMERICAL GRADIENT!")
      }else{
        warning("IT WORKS WITH NUMERICAL GRADIENT!")
      }
    }
    
    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+3)]
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    prec <- coef[p1+2] # precision parameter
    lambda<- coef[p1+3]
    
    z$alpha <- alpha
    z$phi <- phi
    z$prec <- prec
    z$lambda <- lambda
    
    errorhat<-rep(0,n) # E(error)=0 
    etahat<-rep(NA,n)
    
    for(i in (m+1):n)
    {
      etahat[i]<-alpha + (phi%*%ynew[i-ar]) 
      errorhat[i]<-ynew[i]-etahat[i] # predictor scale
      #errorhat[i] <- y[i]-linkinv(etahat[i]) # original scale
      #mui<-linkinv(etahat[i])
      #errorhat[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
    }
    muhat <- linkinv(etahat[(m+1):n],lambda)
    y1<-y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    z$mustarhat <- digamma(muhat * prec) - digamma((1 - muhat) * prec)
    
    #z$resid1 <- errorhat
    
    vI <- as.vector(rep(1,n-m)) 
    
    psi1 = trigamma(muhat * prec)
    psi2 = trigamma((1 - muhat) * prec)
    mT <- diag(mu.eta(etahat[(m+1):n],lambda)) #muhat por etahat[(m+1):n]
    # W = diag(c(prec * (psi1 + psi2))) %*% mT^2
    VT2 = diag(c(prec^2 * (psi1 + psi2))) %*% mT^2
    VT = diag(c(prec^2 * (psi1 + psi2))) %*% mT   #v_t * T
    vc = prec * (psi1 * muhat - psi2 * (1 - muhat))
    D = diag(c(psi1 * (muhat^2) + psi2 * (1 - muhat)^2 - trigamma(prec)))
    rL <- (((1+lambda*exp(etahat[(m+1):n]))^(-1/lambda))/lambda)*( (1/(exp(-etahat[(m+1):n])+lambda))*(1+0)-(log(1+lambda*exp(etahat[(m+1):n]))/lambda) ) 
    V = diag(c(prec^2 * (psi1 + psi2)))   #v_t

    Kaa <- as.matrix(sum(diag(VT2)))
    Kpa <- t(P) %*% VT2 %*% vI
    Kap <- t(Kpa)
    Kaprec <- vI %*% mT %*% vc
    Kpreca <- Kaprec
    Kpp <- t(P) %*% VT2 %*% P
    Kpprec <- t(P) %*% mT %*% vc
    Kprecp <- t(Kpprec)
    Kprecprec <- sum(diag(D))
    
    Kal <- vI %*% VT %*% rL
    Kla <- t(Kal)
    Kpl <- t(P) %*% VT %*% rL
    Klp <- t(Kpl)
    Kprecl <- vc %*% rL
    Klprec <- t(Kprecl)
    Kll <- t(rL) %*% V %*% rL
    
    K <- rbind(
      cbind(Kaa,Kap,Kaprec,Kal),
      cbind(Kpa,Kpp,Kpprec,Kpl),
      cbind(Kpreca,Kprecp,Kprecprec,Kprecl),
      cbind(Kla,Klp,Klprec,Kll)
    )
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha + (phi%*%ynew_prev[n+i-ar]) 
      #errorhat[i]<-ynew[i]-etahat[i] # predictor scale
      y_prev[n+i] <- linkinv(ynew_prev[n+i],lambda)
    }
    
  }
  
  ######### only BMA model
  if(any(is.na(ar)==T) && any(is.na(ma)==F) && any(is.na(X)==T))
  { 
    p1<-0
    print("BMA model",quote=F)
    reg <- c(mqo,rep(0,q1), prec,lambda) # initializing the parameter values
    
    loglik <- function(z) 
    {
      alpha <- z[1]
      theta <- z[2:(q1+1)]
      prec <- z[q1+2] # precision parameter
      lambda <- z[q1+3]
      
      eta <- error <- rep(0,n) # E(error)=0 
      
      for(i in (m+1):n)
      {
        eta[i] <- alpha + (theta%*%error[i-ma])
        error[i]<-ynew[i]-eta[i]
        #error[i] <- y[i]-linkinv(eta[i])
        #mui<-linkinv(eta[i])
        #error[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
      }
      mu <- linkinv(eta[(m+1):n],lambda)
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings(dbeta(y1, mu * prec, (1 - mu) * prec, log = TRUE))
      sum(ll)
    } 
    
    escore <- function(z) 
    {
      alpha <- z[1]
      theta <- z[2:(q1+1)]
      prec <- z[q1+2] # precision parameter
      lambda <- z[q1+3]
      
      error<- rep(0,n) # E(error)=0 
      eta<- rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i] <- alpha + (theta%*%error[i-ma])
        error[i]<- ynew[i]-eta[i] 
        #error[i] <- y[i]-linkinv(eta[i])
        #mui<-linkinv(eta[i])
        #error[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
      }
      
      mu <- linkinv(eta[(m+1):n],lambda)
      y1<-y[(m+1):n]
      ystar <- log(y1/(1-y1))
      mustar <- digamma(mu * prec) - digamma((1 - mu) * prec)
      
      mT <- diag(mu.eta(eta[(m+1):n],lambda)) #mu changed by aquilo
      
      R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
      
      for(i in 1:(n-m))
      {
        R[i,] <- error[i+m-ma]
      }
      
      ###FB  recorrencias
      deta.dalpha <- rep(0,n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)
      deta.dlambda<- rep(0,n) #matrix(0, ncol=1, nrow=n) #me
      dmu.dlambda<- rep(0,n)
      
      # print(dim(P))
      # print(dim(deta.dphi))
      
      for(i in (m+1):n)
      {
        #print(P[(i-m),(i-m)])
        deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
        deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
        
        dmu.dlambda[i] <- (((1+lambda*exp(eta[i]))^(-1/lambda))/lambda)*( (1/(exp(-eta[i])+lambda))*(1+deta.dlambda[i])-(log(1+lambda*exp(eta[i]))/lambda) ) 
      }
      
      a <- deta.dalpha[(m+1):n]
      rR <- deta.dtheta[(m+1):n,]        
      rL <- dmu.dlambda[(m+1):n]

      
      Ualpha <- prec * a %*% mT %*% (ystar-mustar)
      Utheta <- prec * t(rR) %*% mT %*% (ystar-mustar)
      Uprec <-  sum(mu * (ystar-mustar) + log(1-y1) 
                    - digamma((1 - mu) * prec) + digamma(prec) )
      Ulambda <- prec *sum((ystar-mustar) * rL)
      
      rval <- c(Ualpha,Utheta,Uprec,Ulambda)
    }
    names_par <- c("alpha",names_theta,"precision","lambda")
    
    opt <- optim(reg, loglik, escore, method = "BFGS", 
                 control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
    
    if (opt$conv != 0) ###
    {
      warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
      opt <- optim(reg, loglik, #escore,
                   method = "BFGS", 
                   control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
      if (opt$conv != 0)
      {
        warning("FUNCTION DID NOT CONVERGE NEITHER WITH NUMERICAL GRADIENT!")
      }else{
        warning("IT WORKS WITH NUMERICAL GRADIENT!")
      }
    }
    
    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(q1+3)]
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <-coef[1]
    theta <- coef[2:(q1+1)]
    prec <- coef[q1+2] # precision parameter
    lambda <- coef[q1+3]
    
    z$alpha <- alpha
    z$theta <- theta
    z$prec <- prec
    z$lambda <- lambda
    
    errorhat<-rep(0,n) # E(error)=0 
    etahat<-rep(NA,n)
    
    for(i in (m+1):n)
    {
      etahat[i]<-alpha + (theta%*%errorhat[i-ma])
      errorhat[i]<-ynew[i]-etahat[i] # predictor scale
      #errorhat[i] <- y[i]-linkinv(etahat[i]) # original scale
      #mui<-linkinv(etahat[i])
      #errorhat[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
    }
    muhat <- linkinv(etahat[(m+1):n],lambda)
    y1<-y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    z$mustarhat <- digamma(muhat * prec) - digamma((1 - muhat) * prec)
    #z$resid1 <- errorhat
    
    R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
    for(i in 1:(n-m))
    {
      R[i,] <- errorhat[i+m-ma]
    }
    
    vI <- as.vector(rep(1,n-m)) 
    
    ###FB  recorrencias
    deta.dalpha <- rep(0,n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)
    deta.dlambda<- rep(0,n) #matrix(0, ncol=1, nrow=n) #me
    dmu.dlambda<- rep(0,n)
    
    # print(dim(P))
    # print(dim(deta.dphi))
    
    for(i in (m+1):n)
    {
      #print(P[(i-m),(i-m)])
      deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
      deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
      
      dmu.dlambda[i] <- (((1+lambda*exp(etahat[i]))^(-1/lambda))/lambda)*( (1/(exp(-etahat[i])+lambda))*(1+deta.dlambda[i])-(log(1+lambda*exp(etahat[i]))/lambda) )
    }
    
    a <- deta.dalpha[(m+1):n]
    rR <- deta.dtheta[(m+1):n,]        
    rL <- dmu.dlambda[(m+1):n]
    
    psi1 = trigamma(muhat * prec)
    psi2 = trigamma((1 - muhat) * prec)
    mT <- diag(mu.eta(etahat[(m+1):n],lambda)) #muhat por etahat[(m+1):n]
    # W = diag(c(prec * (psi1 + psi2))) %*% mT^2
    VT2 = diag(c(prec^2 * (psi1 + psi2))) %*% mT^2   #v_t * T^2
    VT = diag(c(prec^2 * (psi1 + psi2))) %*% mT   #v_t * T
    V = diag(c(prec^2 * (psi1 + psi2)))   #v_t
    vc = prec * (psi1 * muhat - psi2 * (1 - muhat))
    D = diag(c(psi1 * (muhat^2) + psi2 * (1 - muhat)^2 - trigamma(prec)))
    
    Kaa <- t(a) %*% VT2 %*% a # r ok
    Kta <- t(rR) %*% VT2 %*% a # ok
    Kat <- t(Kta) # ok 
    Kaprec <- t(a) %*% mT %*% vc # ok
    Kpreca <- Kaprec #ok 
    Ktt <- t(rR) %*% VT2 %*% rR
    Ktprec <- t(rR) %*% mT %*% vc
    Kprect <- t(Ktprec)
    Kprecprec <- sum(diag(D))
    
    Kal <- t(a) %*% VT %*% rL
    Kla <- t(Kal)
    Ktl <- t(rR) %*% VT %*% rL #rP nu fodi
    Klt <- t(Ktl)
    Kprecl <- vc %*% rL
    Klprec <- t(Kprecl)
    Kll <- t(rL) %*% V %*% rL
    
    K <- rbind(
      cbind(Kaa,Kat,Kaprec,Kal),
      cbind(Kta,Ktt,Ktprec,Ktl),
      cbind(Kpreca,Kprect,Kprecprec,Kprecl),
      cbind(Kla,Klt,Klprec,Kll)
    )
    
    # print(K)
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha + (theta%*%errorhat[n+i-ma])
      #errorhat[i]<-ynew[i]-etahat[i] # predictor scale
      y_prev[n+i] <- linkinv(ynew_prev[n+i],lambda)
      errorhat[n+i] <- 0 # original scale
    }
  }
  
  
  ##### BAR X model
  if(any(is.na(ar)==F) && any(is.na(ma)==T) && any(is.na(X)==F))
  {
    print("BAR~x model",quote=F)
    q1<-0
    beta1<- mqo[(p1+2):length(mqo)]
    reg <- c(mqo[1:(p1+1)], prec, lambda, beta1) # initializing the parameter values
    
    loglik <- function(z) 
    {
      alpha <- z[1]
      phi = z[2:(p1+1)] 
      prec <- z[p1+2] # precision parameter
      lambda<-z[p1+3]
      beta <- z[(p1+4):length(z)]
      
      error<-rep(0,n) # E(error)=0 
      eta<-rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) 
        error[i] <- ynew[i]-eta[i] 
      }
      mu <- linkinv(eta[(m+1):n],lambda)
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings(dbeta(y1, mu * prec, (1 - mu) * prec, log = TRUE))
      sum(ll)
    } 
    
    escore <- function(z) 
    {
      alpha <- z[1]
      phi <- z[2:(p1+1)]
      prec <- z[p1+2] # precision parameter
      lambda<-z[p1+3]
      beta <- z[(p1+4):length(z)]
      
      error<-rep(0,n) 
      eta<- rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) ))
        error[i] <- ynew[i]-eta[i]       
      }
      
      mu <- linkinv(eta[(m+1):n],lambda)
      y1<-y[(m+1):n]
      ystar <- log(y1/(1-y1))
      mustar <- digamma(mu * prec) - digamma((1 - mu) * prec)
      
      mT <- diag(mu.eta(eta[(m+1):n],lambda)) #changed
      rL <- (((1+lambda*exp(eta[(m+1):n]))^(-1/lambda))/lambda)*( (1/(exp(-eta[(m+1):n])+lambda))*(1+0)-(log(1+lambda*exp(eta[(m+1):n]))/lambda) ) 
      
      P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
      for(i in 1:(n-m))
      {
        P[i,] <- ynew[i+m-ar] - X[i+m-ar,]%*%as.matrix(beta)
      }
      
      M <- matrix(rep(NA,(n-m)*length(beta)),ncol=length(beta))
      for(i in 1:(n-m))
      {
        for(j in 1:length(beta))
          M[i,j] <- X[i+m,j]-sum(phi*X[i+m-ar,j])  
      }
      
      
      Ualpha <- prec * sum((mu.eta(eta[(m+1):n],lambda)) * (ystar-mustar))  ##mu por eta[(m+1):n] sera?
      Uphi <-   prec * t(P) %*% mT %*% (ystar-mustar)
      Uprec <-  sum(mu * (ystar-mustar) + log(1-y1) 
                    - digamma((1 - mu) * prec) + digamma(prec) )
      Ulambda <- prec * sum((ystar-mustar) * rL)
      Ubeta <-   prec * t(M) %*% mT %*% (ystar-mustar)
      
      
      rval <- c(Ualpha,Uphi,Uprec,Ulambda,Ubeta)
    }
    names_par <- c("alpha",names_phi,"precision","lambda",names_beta)
    
    opt <- optim(reg, loglik, escore, method = "BFGS", 
                 control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
    
    if (opt$conv != 0) ###
    {
      warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
      opt <- optim(reg, loglik, #escore,
                   method = "BFGS", 
                   control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
      if (opt$conv != 0)
      {
        warning("FUNCTION DID NOT CONVERGE NEITHER WITH NUMERICAL GRADIENT!")
      }else{
        warning("IT WORKS WITH NUMERICAL GRADIENT!")
      }
    }
    
    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+3+ncol(X))]
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    prec <- coef[p1+2] # precision parameter
    lambda<- coef[p1+3]
    beta <- coef[(p1+4):length(coef)]
    
    
    z$alpha <- alpha
    z$phi <- phi
    z$prec <- prec
    z$lambda <- lambda
    
    errorhat<-rep(0,n) # E(error)=0 
    etahat<-rep(NA,n)
    
    for(i in (m+1):n)
    {
      etahat[i]<-alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) 
      errorhat[i] <- ynew[i]-etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n],lambda)
    y1<-y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    z$mustarhat <- digamma(muhat * prec) - digamma((1 - muhat) * prec)
    
    P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
    for(i in 1:(n-m))
    {
      P[i,] <- ynew[i+m-ar] - X[i+m-ar,]%*%as.matrix(beta)
    }
    
    M <- matrix(rep(NA,(n-m)*length(beta)),ncol=length(beta))
    for(i in 1:(n-m))
    {
      for(j in 1:length(beta))
        M[i,j] <- X[i+m,j]-sum(phi*X[i+m-ar,j])  
    }
    
    vI <- as.vector(rep(1,n-m)) 
    
    psi1 = trigamma(muhat * prec)
    psi2 = trigamma((1 - muhat) * prec)
    mT <- diag(mu.eta(etahat[(m+1):n],lambda)) #muhat por etahat[(m+1):n]
    # W = diag(c(prec * (psi1 + psi2))) %*% mT^2
    VT2 = diag(c(prec^2 * (psi1 + psi2))) %*% mT^2
    VT = diag(c(prec^2 * (psi1 + psi2))) %*% mT   #v_t * T
    vc = prec * (psi1 * muhat - psi2 * (1 - muhat))
    D = diag(c(psi1 * (muhat^2) + psi2 * (1 - muhat)^2 - trigamma(prec)))
    rL <- (((1+lambda*exp(etahat[(m+1):n]))^(-1/lambda))/lambda)*( (1/(exp(-etahat[(m+1):n])+lambda))*(1+0)-(log(1+lambda*exp(etahat[(m+1):n]))/lambda) ) 
    V = diag(c(prec^2 * (psi1 + psi2)))   #v_t
    
    Kaa <- as.matrix(sum(diag(VT2)))
    Kpa <- t(P) %*% VT2 %*% vI
    Kap <- t(Kpa)
    Kaprec <- vI %*% mT %*% vc
    Kpreca <- Kaprec
    Kpp <- t(P) %*% VT2 %*% P
    Kpprec <- t(P) %*% mT %*% vc
    Kprecp <- t(Kpprec)
    Kprecprec <- sum(diag(D))
    
    Kba <- t(M)%*% VT2 %*% vI
    Kbb <- t(M) %*% VT2 %*% M
    Kbprec <- t(M) %*% mT %*% vc 
    Kbp <- t(M) %*% VT2 %*% P
    
    Kab <- t(Kba)
    Kprecb <- t(Kbprec)
    Kpb <- t(Kbp)
    
    Kal <- vI %*% VT %*% rL
    Kla <- t(Kal)
    
    Kbl <- t(M) %*% VT %*% rL
    Klb <- t(Kbl)
    
    Kpl <- t(P) %*% VT %*% rL
    Klp <- t(Kpl)
    Kprecl <- vc %*% rL
    Klprec <- t(Kprecl)
    Kll <- t(rL) %*% V %*% rL
    
    K <- rbind(
      cbind(Kaa,Kap,Kaprec,Kal,Kab),
      cbind(Kpa,Kpp,Kpprec,Kpl,Kpb),
      cbind(Kpreca,Kprecp,Kprecprec,Kprecl,Kprecb),
      cbind(Kla,Klp,Klprec,Kll,Klb),
      cbind(Kba,Kbp,Kbprec,Kbl,Kbb)
    )
    # print(K)
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    X_prev<- rbind(X,X_hat)
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) + (phi%*%(ynew_prev[n+i-ar]-X_prev[n+i-ar,]%*%as.matrix(beta) ))
      #errorhat[i]<-ynew[i]-etahat[i] # predictor scale
      y_prev[n+i] <- linkinv(ynew_prev[n+i],lambda)
      errorhat[n+i] <- 0 # original scale
    }
    
  }
  
  ######### BMAX model
  if(any(is.na(ar)==T) && any(is.na(ma)==F) && any(is.na(X)==F))
  { 
    p1<-0
    print("BMA~X model",quote=F)
    beta1<- mqo[(2):length(mqo)]
    reg <- c(mqo[1],rep(0,q1), prec,lambda,beta1) # initializing the parameter values
    
    loglik <- function(z) 
    {
      alpha <- z[1]
      theta <- z[2:(q1+1)]
      prec <- z[q1+2] # precision parameter
      lambda <- z[q1+3]
      beta <- z[(q1+4):length(z)]
      
      error<-rep(0,n) # E(error)=0 
      eta<-rep(NA,n)      
      for(i in (m+1):n)
      {
        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (theta%*%error[i-ma])
        error[i] <- ynew[i]-eta[i] 
      }
      mu <- linkinv(eta[(m+1):n],lambda)
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings(dbeta(y1, mu * prec, (1 - mu) * prec, log = TRUE))
      sum(ll)
    } 
    
    escore <- function(z) 
    {
      alpha <- z[1]
      theta <- z[2:(q1+1)]
      prec <- z[q1+2] # precision parameter
      lambda <- z[q1+3]      
      beta <- z[(q1+4):length(z)]

      
      error<- rep(0,n) # E(error)=0 
      eta<- rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i]<- alpha + X[i,]%*%as.matrix(beta) + (theta%*%error[i-ma])
        error[i] <- ynew[i]-eta[i] 
      }
      
      mu <- linkinv(eta[(m+1):n],lambda)
      y1<-y[(m+1):n]
      ystar <- log(y1/(1-y1))
      mustar <- digamma(mu * prec) - digamma((1 - mu) * prec)
      
      mT <- diag(mu.eta(eta[(m+1):n],lambda)) #mu changed by aquilo
      
      R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
      for(i in 1:(n-m))
      {
        R[i,] <- error[i+m-ma]
      }
      
      k1<- length(beta)
      M <- matrix(rep(NA,(n-m)*k1),ncol=k1)
      for(i in 1:(n-m))
      {
        for(j in 1:k1)
          M[i,j] <- X[i+m,j] 
      }
      
      ###FB  recorrencias
      deta.dalpha <- rep(0,n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)
      deta.dlambda<- rep(0,n) #matrix(0, ncol=1, nrow=n) #me
      dmu.dlambda<- rep(0,n)
      deta.dbeta<- matrix(0, ncol=k1,nrow=n)
      
      for(i in (m+1):n)
      {
        #print(P[(i-m),(i-m)])
        deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
        deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
        deta.dbeta[i,]<- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
        
        dmu.dlambda[i] <- (((1+lambda*exp(eta[i]))^(-1/lambda))/lambda)*( (1/(exp(-eta[i])+lambda))*(1+deta.dlambda[i])-(log(1+lambda*exp(eta[i]))/lambda) ) 
      }
      
      a <- deta.dalpha[(m+1):n]
      rR <- deta.dtheta[(m+1):n,]        
      rL <- dmu.dlambda[(m+1):n]
      rM <- deta.dbeta[(m+1):n,]
      
      Ualpha <- prec * a %*% mT %*% (ystar-mustar)
      Utheta <- prec * t(rR) %*% mT %*% (ystar-mustar)
      Uprec <-  sum(mu * (ystar-mustar) + log(1-y1) 
                    - digamma((1 - mu) * prec) + digamma(prec) )
      Ulambda <- prec *sum((ystar-mustar) * rL)
      Ubeta <-   prec * t(rM) %*% mT %*% (ystar-mustar)
      
      rval <- c(Ualpha,Utheta,Uprec,Ulambda,Ubeta)
    }
    names_par <- c("alpha",names_theta,"precision","lambda",names_beta)
    
    opt <- optim(reg, loglik, escore, method = "BFGS", 
                 control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
    
    if (opt$conv != 0) ###
    {
      warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
      opt <- optim(reg, loglik, #escore,
                   method = "BFGS", 
                   control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
      if (opt$conv != 0)
      {
        warning("FUNCTION DID NOT CONVERGE NEITHER WITH NUMERICAL GRADIENT!")
      }else{
        warning("IT WORKS WITH NUMERICAL GRADIENT!")
      }
    }
    
    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(q1+3+ncol(X))]
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <-coef[1]
    theta <- coef[2:(q1+1)]
    prec <- coef[q1+2] # precision parameter
    lambda <- coef[q1+3]
    beta <- coef[(q1+4):length(coef)]
    
    z$alpha <- alpha
    z$theta <- theta
    z$prec <- prec
    z$lambda <- lambda
    
    errorhat<-rep(0,n) # E(error)=0 
    etahat<-rep(NA,n)
    
    for(i in (m+1):n)
    {
      etahat[i]<-alpha + X[i,]%*%as.matrix(beta) + (theta%*%errorhat[i-ma])
      errorhat[i] <- ynew[i]-etahat[i] # 
    }
    muhat <- linkinv(etahat[(m+1):n],lambda)
    y1<-y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    z$mustarhat <- digamma(muhat * prec) - digamma((1 - muhat) * prec)
    #z$resid1 <- errorhat
    
    R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
    for(i in 1:(n-m))
    {
      R[i,] <- errorhat[i+m-ma]
    }
    
    k1<-length(beta)
    M <- matrix(rep(NA,(n-m)*k1),ncol=k1)
    for(i in 1:(n-m))
    {
      for(j in 1:k1)
        M[i,j] <- X[i+m,j]
    }
    ###FB  recorrencias
    deta.dalpha <- rep(0,n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)
    deta.dlambda<- rep(0,n) #matrix(0, ncol=1, nrow=n) #me
    dmu.dlambda<- rep(0,n)
    deta.dbeta<- matrix(0, ncol=k1,nrow=n)

    for(i in (m+1):n)
    {
      #print(P[(i-m),(i-m)])
      deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
      deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
      
      dmu.dlambda[i] <- (((1+lambda*exp(etahat[i]))^(-1/lambda))/lambda)*( (1/(exp(-etahat[i])+lambda))*(1+deta.dlambda[i])-(log(1+lambda*exp(etahat[i]))/lambda) )
      deta.dbeta[i,]<- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
      
      }
    
    a <- deta.dalpha[(m+1):n]
    rR <- deta.dtheta[(m+1):n,]        
    rL <- dmu.dlambda[(m+1):n]
    rM <- deta.dbeta[(m+1):n,]
    
    
    psi1 = trigamma(muhat * prec)
    psi2 = trigamma((1 - muhat) * prec)
    mT <- diag(mu.eta(etahat[(m+1):n],lambda)) #muhat por etahat[(m+1):n]
    # W = diag(c(prec * (psi1 + psi2))) %*% mT^2
    VT2 = diag(c(prec^2 * (psi1 + psi2))) %*% mT^2   #v_t * T^2
    VT = diag(c(prec^2 * (psi1 + psi2))) %*% mT   #v_t * T
    V = diag(c(prec^2 * (psi1 + psi2)))   #v_t
    vc = prec * (psi1 * muhat - psi2 * (1 - muhat))
    D = diag(c(psi1 * (muhat^2) + psi2 * (1 - muhat)^2 - trigamma(prec)))
    
    Kaa <- t(a) %*% VT2 %*% a # r ok
    Kta <- t(rR) %*% VT2 %*% a # ok
    Kat <- t(Kta) # ok 
    Kaprec <- t(a) %*% mT %*% vc # ok
    Kpreca <- Kaprec #ok 
    Ktt <- t(rR) %*% VT2 %*% rR
    Ktprec <- t(rR) %*% mT %*% vc
    Kprect <- t(Ktprec)
    Kprecprec <- sum(diag(D))
    
    Kba <- t(rM)%*% VT2 %*% a # t(a) %*% VT2 %*% rM
    Kbb <-  t(rM) %*% VT2 %*% rM
    Kbprec <- t(rM) %*% mT %*% vc
    Kbt <- t(rM) %*% VT2 %*% rR
    
    Kab <- t(Kba)
    Kprecb <- t(Kbprec)
    Ktb <- t(Kbt)
    
    Kal <- t(a) %*% VT %*% rL
    Kla <- t(Kal)
    Kbl <- t(rM) %*% VT %*% rL
    Klb <- t(Kbl)
    Ktl <- t(rR) %*% VT %*% rL #rP nu fodi
    Klt <- t(Ktl)
    Kprecl <- vc %*% rL
    Klprec <- t(Kprecl)
    Kll <- t(rL) %*% V %*% rL
    
    K <- rbind(
      cbind(Kaa,Kat,Kaprec,Kal,Kab),
      cbind(Kta,Ktt,Ktprec,Ktl,Ktb),
      cbind(Kpreca,Kprect,Kprecprec,Kprecl,Kprecb),
      cbind(Kla,Klt,Klprec,Kll,Klb),
      cbind(Kba,Kbt,Kbprec,Kbl,Kbb)
    )
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    X_prev<- rbind(X,X_hat)#X[i - ar, ]
    
    #     (phi%*%(ynew[i-ar]-X[i-ar,]%*%beta ))
    #     (phi%*%(ynew_prev[n+i-ar]-X_prev[i-ar,]%*%beta ))
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) + (theta%*%errorhat[n+i-ma])
      #errorhat[i]<-ynew[i]-etahat[i] # predictor scale
      y_prev[n+i] <- linkinv(ynew_prev[n+i],lambda)
      errorhat[n+i] <- 0 # original scale
    }
  }
  
  ##############################################
  { 
    # residuals
    res1 <- y-z$fitted #padronizado
    z$residP <- (res1[(m+1):n])/sqrt(z$fitted[(m+1):n]*(1-z$fitted[(m+1):n])/(1+z$prec))
    
    # res2 <- ynew-z$etahat
    # z$resid2 <- (res2[(m+1):n])/sqrt((mu.eta(z$fitted[(m+1):n])^{-2})*(z$fitted[(m+1):n]*(1-z$fitted[(m+1):n])/(1+z$prec)) )
    # com lambda z$resid2 <- (res2[(m+1):n])/sqrt((mu.eta(z$fitted[(m+1):n],lambda)^{-2})*(z$fitted[(m+1):n]*(1-z$fitted[(m+1):n])/(1+z$prec)) )
   
    res2 <- as.vector(ystar)[(m+1):n]-z$mustarhat # ponderado
    z$residW <- (res2/sqrt( trigamma(z$fitted[(m+1):n]*z$prec)+trigamma((1-z$fitted[(m+1):n])*z$prec) ) )

    l_tilde <- (dbeta(y[(m+1):n], y[(m+1):n] * z$prec, (1 - y[(m+1):n]) * z$prec, log = TRUE))
    l_hat <- (dbeta(y[(m+1):n], z$fitted[(m+1):n] * z$prec, (1 - z$fitted[(m+1):n]) * z$prec, log = TRUE))

    dt <- (l_tilde-l_hat)
    dt[which(dt<0)]<-0

    z$residD <- sign(y[(m+1):n]-z$fitted[(m+1):n])*sqrt(2*(dt)) #deviance
    
    z$residQ <- qnorm(pbeta(y[(m+1):n], z$fitted[(m+1):n]*z$prec, (1-z$fitted[(m+1):n])*z$prec))
    
    
    if(resid==1) residc <- z$residP
    if(resid==2) residc <- z$residW
    if(resid==3) residc <- z$residD
    if(resid==4) residc <- z$residQ

    # print(Box.test(z$residP, lag = 20, type =  "Ljung-Box", fitdf = p1+q1+3+length(beta)))
    # print(Box.test(z$residW, lag = 20, type =  "Ljung-Box", fitdf = p1+q1+3+length(beta)))
    # print(Box.test(z$residD, lag = 20, type =  "Ljung-Box", fitdf = p1+q1+3+length(beta)))
    # print(Box.test(z$residQ, lag = 20, type =  "Ljung-Box", fitdf = p1+q1+3+length(beta)))
    

    ############################################

    #vcov <- solve(K)
    vcov <- chol2inv(chol(K))
    z$vcov <- vcov

    stderror <- sqrt(diag(vcov)) 
    z$stderror <- stderror 

    z$zstat <- abs(z$coef/stderror)
    z$pvalues <- 2*(1 - pnorm(z$zstat))

    #z$loglik <- opt$value
    z$loglik <- opt$value*(n/(n-m))
    z$counts <- as.numeric(opt$counts[1])
    #   z$aic <- -2*z$loglik+2*(p1+q1+2)
    #   z$bic <- -2*z$loglik+log(n)*(p1+q1+2)
    
    if(any(is.na(X)==F))
    {
      z$aic <- -2*z$loglik+2*(p1+q1+3+length(beta))#+3
      z$bic <- -2*z$loglik+log(n)*(p1+q1+3+length(beta))#+3
    }else{
      z$aic <- -2*z$loglik+2*(p1+q1+2)
      z$bic <- -2*z$loglik+log(n)*(p1+q1+2)
    }
  } 
  
  
  z$serie <- y
  z$barma <- names_par
  
  
  model_presentation <- cbind(round(z$coef,3),round(z$stderror,3),round(z$zstat,3),round(z$pvalues,3))
  colnames(model_presentation)<-c("Estimate","Std. Error","z value","Pr(>|z|)")
  
  z$model <- model_presentation
  # z$link <- link
  
  z$forecast <- y_prev[(n+1):(n+h1)]
  z$y_prev <- y_prev

  

  
  ######## GRAPHICS ################################
  {

      print(model_presentation)
      print(" ",quote=F)
  print(c("Log-likelihood:",round(z$loglik,4)),quote=F)
      print(c("Number of iterations in BFGS optim:",z$counts),quote=F)
  print(c("AIC:",round(z$aic,4)," BIC:",round(z$bic,4)),quote=F)

      print("Residuals:",quote=F)
      print(summary(residc))

      t<-seq(-5,n+6,by=1)

      par(mfrow=c(1,1))
      par(mar=c(2.8, 2.7, 1.2, 1)) # margens c(baixo,esq,cima,direia)
      par(mgp=c(1.7, 0.45, 0))
      plot(residc,main=" ",xlab="Index",ylab="Residuals", pch = "+",ylim=c(-4,4))
      lines(t,rep(-3,n+12),lty=2,col=1)
      lines(t,rep(3,n+12),lty=2,col=1)
      lines(t,rep(-2,n+12),lty=3,col=1)
      lines(t,rep(2,n+12),lty=3,col=1)


      max_y<- max(c(z$fitted,y),na.rm=T)
      min_y<- min(c(z$fitted,y),na.rm=T)
      plot(as.vector(z$fitted), as.vector(y), main=" ", pch = "+",
           xlab="Fitted values",ylab="Observed values",
           xlim=c(0.95*min_y,max_y*1.05),
           ylim=c(0.95*min_y,max_y*1.05))
      lines(c(-0.2,1.2),c(-0.2,1.2),lty=2)

      densidade<-density(residc)
      plot(densidade,ylab=" ",xlab=" ",main=" ")
      lines(densidade$x,dnorm(densidade$x),lty=2)
      legend("topright",c("Estimated density","Standard normal"),#pch=vpch,
             pt.bg="white", lty=c(1,2), bty="n")

      acf(residc,ylab="FAC",xlab="lag") # função de autocorrelação

      pacf(residc,ylab="FACP",xlab="lag") # função de autocorrelação parcial

      max_r<- max(residc,na.rm=T)
      min_r<- min(residc,na.rm=T)
      qqnorm(residc, pch = "+",
             xlim=c(0.95*min_r,max_r*1.05),
             ylim=c(0.95*min_r,max_r*1.05),
             main="",xlab="Normal quantiles",ylab="Empirical quantiles")
      lines(c(-10,10),c(-10,10),lty=2)

      par(mfrow=c(1,1))
      plot(y,type="l",ylab="Relative air humidity",xlab="Time")
      lines(z$fitted,type = "l",lty=2)
      legend(2007.1,.69, lty = c(1,2), cex=0.96 ,bty ="n", 
             legend = c("Observed values","Fitted values"))
      
             # legend = c(expression(paste(beta,"ARMA(3,3) proposto"),paste(beta,"ARMA(3,3) usual"))))
      
      fim<-end(y)[1]+end(y)[2]/12

      y_prev <- ts(y_prev, start=start(y), frequency=frequency(y))
      par(mfrow=c(1,1))
      plot(y_prev,type="l",col="red", ylim=c(min(y),max(y)),ylab="Relative humidity",xlab="Time")
      abline(v=fim,lty=2)
      lines(y)

      w1<-3
      h1<-3
      
      if(diag>1)
      { 

        pdf(file = "resid_v_ind.pdf",width = w1, height = h1,family = "Times")
        {
          par(mfrow=c(1,1))
          par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
          par(mgp=c(1.7, 0.45, 0))
          plot(residc,main=" ",xlab="Indices",ylab="residuos", pch = "+",ylim=c(-4,4))
          lines(t,rep(-3,n+12),lty=2,col=1)
          lines(t,rep(3,n+12),lty=2,col=1)
          lines(t,rep(-2,n+12),lty=3,col=1)
          lines(t,rep(2,n+12),lty=3,col=1)
        }
        dev.off()

        pdf(file = "obs_v_fit.pdf",width = w1, height = h1,family = "Times")
        {
          par(mfrow=c(1,1))
          par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
          par(mgp=c(1.7, 0.45, 0))
          plot(as.vector(z$fitted), as.vector(y), main=" ", pch = "+",
               xlab="Fitted values",ylab="Observed values",
               xlim=c(0.95*min_y,max_y*1.05),
               ylim=c(0.95*min_y,max_y*1.05))
          lines(c(-0.2,1.2),c(-0.2,1.2),lty=2)
        }
        dev.off()

        pdf(file = "resid_density.pdf",width = w1, height = h1,family = "Times")
        {
          par(mfrow=c(1,1))
          par(mar=c(1.5, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
          par(mgp=c(1.7, 0.45, 0))

          plot(densidade,ylab="densidade",main=" ",xlab=" ",ylim=c(0,1.15*max(densidade$y)))
          lines(densidade$x,dnorm(densidade$x),lty=2)
          legend("topleft",c("Estimated density","Standard normal"),#pch=vpch,
                 pt.bg="white", lty=c(1,2), bty="n")
        }
        dev.off()

        pdf(file = "resid_FAC.pdf",width = w1, height = h1,family = "Times")
        {
          par(mfrow=c(1,1))
          par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
          par(mgp=c(1.7, 0.45, 0))
          acf(residc,ylab="FAC",xlab="lag") # função de autocorrelação
        }
        dev.off()

        pdf(file = "resid_FACP.pdf",width = w1, height = h1,family = "Times")
        {
          par(mfrow=c(1,1))
          par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
          par(mgp=c(1.7, 0.45, 0))
          pacf(residc,ylab="FACP",xlab="lag") # função de autocorrelação parcial
        }
        dev.off()

        pdf(file = "qq_plot.pdf",width = w1, height = h1,family = "Times")
        {
          par(mfrow=c(1,1))
          par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
          par(mgp=c(1.7, 0.45, 0))
          qqnorm(residc, pch = "+",
                 xlim=c(0.95*min_r,max_r*1.05),
                 ylim=c(0.95*min_r,max_r*1.05),
                 main="",xlab="Normal quantiles",ylab="Empirical quantiles")
          lines(c(-10,10),c(-10,10),lty=2)
        }
        dev.off()

        pdf(file = "adjusted.pdf",width = 6, height = 4,family = "Times")
        {
          par(mfrow=c(1,1))
          par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
          par(mgp=c(1.7, 0.45, 0))
          plot(y,type="l",ylab="Relative air humidity",xlab="Time")
          lines(z$fitted,col="red")
        }
        dev.off()
        
        pdf(file = "forecast.pdf",width = 6, height = 4,family = "Times")
        {
          par(mfrow=c(1,1))
          par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
          par(mgp=c(1.7, 0.45, 0))
          plot(y_prev,type="l",col="red", ylim=c(min(y),max(y)),ylab="Relative air humidity",xlab="Time")
          abline(v=fim,lty=2)
          lines(y)
        }
        dev.off()
      }
  }
  
  return(z)
}

