
rm(list = ls())
# SIMULAMOS UNA TRIANGULAR ------------------------------------------------


set.seed(pi)
n <- 1000
A <- runif(n,-pi/2,pi/2)
B <- runif(n,-pi/2,pi/2)
X <- A+B
plot(density(X))

f <- Vectorize(function(x) (1/pi^2*x+1/pi)*(x>-pi&&x<0)+(-1/pi^2*x+1/pi)*(x>0&&x<pi) )

f <- Vectorize(\(x) as.numeric((x>=-.5)*(x<=.5))*.5)


curve(f(x),add=F,col=2, lwd = 2, n = 1000, -pi, pi)



# ------------------------------------------------------------------------

curve(f, -4, 4, col =2, lwd= 2, n = 200)


fourier <- function(funcion,K) {
  a <- sapply(0:K, function(k) integrate(function(x) do.call(funcion,list(x))*cos(k*x)/pi,-pi,pi)$value)
  b <- sapply(0:K, function(k) integrate(function(x) do.call(funcion,list(x))*sin(k*x)/pi,-pi,pi)$value)
  a[1] <- a[1]/2
  f <- Vectorize(function(x) sum(a*cos(0:K*x)+b*sin(0:K*x)))
  out <- list(f=f,a=a,b=b)
  class(out) <- 'fou'
  return(out)
}

fourier(f,K = 3)

curve(f(x),-pi,pi,n=500, col  = 2, lwd = 2)
curve(fourier(f,4)$f(x),add=T,col=4, lwd = 2)


library(ggplot2)
library(gridExtra)
theme_set(theme_bw())

p <- ggplot() + 
  xlim(c(-1,1)) + 
  geom_function(fun = f, n = 301, lwd = .8)

grafs <- sapply(c(1,3,5,8:10,20,30,40),function(k){
  p + geom_function(fun = fourier(f,k)$f, aes(col = "fourier"), lwd = .9) + 
    ggtitle(sprintf("Serie truncada con %s coeficientes", k)) +
    theme(legend.position = "none", plot.title = element_text(size=10))
})

do.call("grid.arrange", args = grafs)

par(mfrow=c(3,3),mai=c(.4,.3,.3,.1))
invisible(sapply(c(1:6,c(10,15,20)),function(k){
  curve(fourier(f,k)$f(x),col=2, lwd = 2, -1, 1)
  curve(f(x),-pi,pi,n=500,main=paste0('k=',k), add = T)
  # curve(f(x),-pi,pi,n=500,main=paste0('k=',k))
  # curve(fourier(f,k)$f(x),add=T,col=2, lwd = 2)
}))



plot.fou <- function(x, ...){plot(x$f,...)}
print.fou <- function(x, ...){print(list(a=round(x$a,3),b=round(x$b,3)),...)}

dev.off()

out <- fourier(f,10)
plot(f,-pi,pi)
plot(out,-pi,pi,add=T,col=2)
curve(out$f(x),-pi,pi,add=T,col=2)


#############

fourier_est <- function(X,K) {
  a <- sapply(0:K, function(k) mean(cos(k*X))/pi)
  b <- sapply(0:K, function(k) mean(sin(k*X))/pi)
  a[1] <- a[1]/2
  f <- Vectorize(function(x) sum(a*cos(0:K*x)+b*sin(0:K*x)))
  out <- list(f=f,a=a,b=b)
  class(out) <- 'fou'
  return(out)
}


mean(cos(2*X)/pi)



# vamos -------------------------------------------------------------------
# set.seed(pi)

n <- 1e4
A <- runif(n,-pi/2,pi/2)
B <- runif(n,-pi/2,pi/2)
X <- A+B
N <- 10


X <- runif(n, -.5, .5)

plot(f,-pi,pi, ylim = c(0,.34))
curve(fourier(f, N)$f(x),add=T,col=2, lwd = 2)
plot(fourier_est(X, N), add = T, col = 4, -pi, pi, lwd = 2)

fourier(f, N)
fourier_est(X, N)
# -------------------------------------------------------------------------



par(mfrow=c(3,3),mai=c(.6,.3,.3,.1))
invisible(sapply(7:15,function(k){
  n <- 1e5
  set.seed(1234)
  # X <- runif(n,-pi/2,pi/2)+runif(n,-pi/2,pi/2)
  X <- runif(n, -.5, .5)
  algo <- fourier_est(X,k)
  plot(f,-pi,pi,n=500,main=paste0('k=',k))
  plot(algo,-pi,pi,add=T,col=2,n=1000, lwd=2)
}))

par(mfrow = c(1,1))
X <- runif(n,-pi/2,pi/2)+runif(n,-pi/2,pi/2)
X <- runif(n, -.5, .5)

out <- fourier_est(X,500)
# plot(f,-pi,pi)
# plot(out,-pi,pi,add=T,col=2,n=1001)

library(ggplot2)
ggplot() + 
  xlim(c(-pi,pi)) +
  geom_function(fun = out$f, aes(col = "estimada")) + 
  geom_function(fun = f, aes(col = "teórica")) + 
  theme_bw() 

n <- 10000
X <- runif(n,-pi/2,pi/2)+runif(n,-pi/2,pi/2)
out <- fourier_est(X,500)

library(ggplot2)
ggplot() + 
  xlim(c(-pi,pi)) +
  geom_function(fun = out$f, aes(col = "estimada")) + 
  geom_function(fun = f, aes(col = "teórica"), alpha = .5) + 
  geom_density(aes(x = X, col = "KDE")) +
  theme_bw() 

# plot(out,-pi,pi,lwd = 2, col=2,n=1001)
# curve(f,-pi,pi, add = T)
# lines(density(X), col = 4, lwd = 2)


### HACER LA COMPARACION ENTRE NUCLEOS Y ESTO, VARIANDO EL K Y LA CANTIDAD DE DATOS
### PROBAR CAMBIAR UN DATO DE LA MUESTRA



# k vecinos ---------------------------------------------------------------


library(TDA)
library(mvtnorm)

set.seed(pi)
datos <- rmvnorm(10000,sigma=matrix(c(1,.5,.5,1),2))

z <- outer(X=seq(-3,3,,100),
           Y = seq(-3,3,,100),
           FUN = Vectorize(function(x,y) {dmvnorm(x=c(x,y),mean=c(0,0),sigma=matrix(c(1,.5,.5,1),2))}
                           ,vectorize.args = c('x','y')))

plot(datos,pch=20,cex=.5,xlim=c(-3,3),ylim=c(-3,3))
contour(x=seq(-3,3,,100),y=seq(-3,3,,100),z,col=2,add=T,lwd=2)

seq(-3,3,,100)
plot(expand.grid(x=seq(-3,3,,100),y=seq(-3,3,,100)),pch=20)


est <- knnDE(datos,Grid = expand.grid(x=seq(-3,3,,100),y=seq(-3,3,,100)),k=100)
est <- cbind(expand.grid(x=seq(-3,3,,100),y=seq(-3,3,,100)),est)
contour(x=seq(-3,3,,100),y=seq(-3,3,,100),z,col='tomato',levels=seq(0.02,0.18,0.02))
contour(x=seq(-3,3,,100),y=seq(-3,3,,100),z=as.matrix(reshape2::dcast(data=est,formula = x~y))[,-1],add=T,levels=seq(0.02,0.18,0.02))

contour(x=seq(-3,3,,100),y=seq(-3,3,,100),z=as.matrix(reshape2::dcast(data=est,formula = x~y))[,-1],add=T,levels=seq(0.02,0.18,0.02))

library(rgl)
persp3d(x=seq(-3,3,,100),y=seq(-3,3,,100),z=as.matrix(reshape2::dcast(data=est,formula = x~y))[,-1],col='tomato',xlab = '',ylab='',zlab='')




# histograma bivariado ----------------------------------------------------



library(mvtnorm)
set.seed(pi)
n1 <- 500
n2 <- 500
n <- n1+n2
datos <- rbind(rmvnorm(n1,sigma=matrix(c(1,-.5,-.5,1),2)),rmvnorm(n2,mean=c(3,3),sigma=matrix(c(1,-.5,-.5,1),2)))

plot(datos,pch=20)

h <- .5

xmin <- -4; xmax <- 7
bins <- seq(xmin,xmax,h)

x0 <- bins[-length(bins)]
x1 <- bins[-1]

datos[1:6,]
cont <- table(apply(apply(datos,1,cut,breaks=bins,labels=x0),2,paste,collapse='_'))/n
library(data.table)
library(magrittr)
as.data.table(cont,keep.rownames = T) %>% 
  .[,c('x','y'):=tstrsplit(V1,'_',fixed=T)] %>% 
  .[,!'V1'] %>%
  .[,(2:3):=lapply(.SD,as.numeric),.SDcols=2:3] %>%
  dcast(x~y,value.var='N',fill=0) %>% 
  .[,!'x'] %>% 
  as.matrix() %>% 
  image
