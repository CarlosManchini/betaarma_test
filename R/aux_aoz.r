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

# diflink <- function(mu,lambda) {
#   ret<-((lambda*(1-mu)^(-(lambda+1)))/((1-mu)^(-lambda)-1))
#   as.vector(ret)
# }
