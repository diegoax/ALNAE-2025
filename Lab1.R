rm(list=ls())

library(ggplot2)
library(gridExtra)
library(scales)
library(data.table)
library(magrittr)

reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))}


# -------------------------------------------------------------------------


n <- 100
h <- .5

xmin <- -4
xmax <- 4
bins <- seq(xmin,xmax,h)

x0 <- bins[-length(bins)]
x1 <- bins[-1]


ggplot(data=data.frame(x=c(-4,4)),aes(x))+stat_function(fun=dnorm)+
  geom_point(data=data.frame(bins),aes(x=bins,y=bins*0),pch='|',size=2)

simular <- function(n,
                    probs = c(.2, .5, .3),
                    medias = c(-1.5, -.2, 1),
                    desvio = c(.3, .2, .2)) {
  z <- rmultinom(1,n,probs)[,1]
  rnorm(n, rep(medias, z), rep(desvio, z))
}



x <- rnorm(n)
cont <- table(cut(x,breaks = bins))/n
f_est <- stepfun(bins[2:(length(bins)-1)],cont/h)

ggplot()+
  geom_point(data=data.frame(x),aes(x,0,col='z'),pch='|',size=1.5,show.legend=F)+
  stat_function(fun=f_est,data=data.frame(x=c(-4,4)),aes(x,col='estimada'),show.legend=F,n = 501)+
  stat_function(fun=dnorm,data=data.frame(x=c(-4,4)),aes(x,col='teórica'),show.legend=F,col=1)


# -------------------------------------------------------------------------



set.seed(pi)

plots <- lapply(seq_len(9),function(eee) {

  x <- rnorm(n)
  cont <- table(cut(x,breaks = bins))/n
  f_est <- stepfun(bins[2:(length(bins)-1)],cont/h)
  grafico <- ggplot()+
    geom_point(data=data.frame(x),aes(x,0,col='z'),pch='|',size=1.5,show.legend=F)+
    # geom_ribbon(data=data.frame(xx=seq(-4,4,length.out = 501)),aes(x=xx,ymin=f_est(xx),ymax=dnorm(xx)),fill='grey',alpha=.6)+
    stat_function(fun=f_est,data=data.frame(x=c(-4,4)),aes(x,col='estimada'),show.legend=F,n = 501)+
    stat_function(fun=dnorm,data=data.frame(x=c(-4,4)),aes(x,col='teórica'),show.legend=F,col=1) + 
    theme_bw()
  
  grafico$funcion <- f_est
  return(grafico)
})

do.call('grid.arrange',plots)
# marrangeGrob(plots, ncol = 1, nrow = 9)


# -------------------------------------------------------------------------



p <- -.25 

do.call('grid.arrange',
        lapply(plots, function(grafico){
          grafico+
            geom_point(data=data.frame(x=p),aes(x,dnorm(x),col='punto'),size=.7)+
            geom_point(data=data.frame(x=p),aes(x,grafico$funcion(x),col='punto'),size=.7)+
            geom_segment(data=data.frame(x=p),aes(x=x,xend=x,y=dnorm(x),yend=grafico$funcion(x),col='punto'),lty=2)+
            theme(legend.position='none')
}))

sapply(plots, function(x) {
  x$funcion(p)
})


# MISE --------------------------------------------------------------------


do.call('grid.arrange',
        lapply(plots, function(grafico){
          grafico+
            geom_ribbon(data=data.frame(xx=seq(-4,4,length.out = 501)),aes(x=xx,ymin=grafico$funcion(xx),ymax=dnorm(xx)),fill='grey',alpha=.6)+
            geom_point(data=data.frame(x=p),aes(x,dnorm(x),col='punto'),size=.7)+
            geom_point(data=data.frame(x=p),aes(x,grafico$funcion(x),col='punto'),size=.7)+
            geom_segment(data=data.frame(x=p),aes(x=x,xend=x,y=dnorm(x),yend=grafico$funcion(x),col='punto'),lty=2)+
            theme(legend.position='none')
        }))


# -------------------------------------------------------------------------


set.seed(pi)

estimadores <- lapply(seq_len(50), function(eeeeeeeee) {
  
  x <- rnorm(n)
  cont <- table(cut(x,breaks = bins))/n
  f_est <- stepfun(bins[2:(length(bins)-1)],cont/h)
  return(f_est)
  
})

# mean_f <- Vectorize(function(p) mean(sapply(estimadores, function(x) x(p))))
# var_f <- Vectorize(function(p) var(sapply(estimadores, function(x) x(p))))
# curve(mean_f, -4, 4, n = 501)
# curve(var_f, -4, 4, n = 501)


base <- ggplot(data.frame(x=c(-4,4)),aes(x))+
  stat_function(fun=dnorm,col=1)

eh <- invisible(lapply(estimadores,function(est){
  base <<- base+
    stat_function(fun=est,aes(col='a'),n=501,show.legend=F,alpha=.3)
}))

eh[[50]]

do.call('grid.arrange',eh[c(2,10,50)])

# -------------------------------------------------------------------------


bins

pnorm(bins[2])-pnorm(bins[1])
pnorm(bins[9])-pnorm(bins[8])

curve(dnorm,-4,4)
curve(dnorm,-.5,0,type='h',col='tomato',add=T)
curve(pnorm,-4,4)
abline(v=bins[8:9],col='tomato')

# -------------------------------------------------------------------------
library(data.table)
library(magrittr)

t <- seq(xmin,xmax,length.out = 201)

set.seed(pi)
estimaciones <- setDT(lapply(seq_len(500),function(foo) {
  x <- rnorm(n)
  cont <- table(cut(x,breaks = bins))/n
  f_est <- stepfun(bins[2:(length(bins)-1)],cont/h)
  f_est(t)
}))

cbind(t,estimaciones) %>% 
  melt(id.vars='t') %>% 
  ggplot(aes(x=t,y=value,group=variable,col=variable))+
    geom_step(show.legend = F,alpha=.05)+
    stat_function(aes(t),fun=dnorm,col=1)


(p1 <- cbind(t,estimaciones) %>% 
  melt(id.vars='t') %>% 
  ggplot()+
  # geom_bin2d(aes(x=t,y=value,group=t,fill=..density..),show.legend=F,binwidth=c(h,1/n/h))+
    geom_bin2d(aes(x=t,y=value,group=t,fill=..density..),show.legend=F,binwidth=c(.04,.04))+
  scale_fill_gradient(name='Densidad',trans=reverselog_trans(10))+
  stat_function(fun=dnorm,show.legend=F) + theme_bw())

# -------------------------------------------------------------------------

  
cbind(t,estimaciones) %>% 
  melt(id.vars='t') %>% 
  .[,bin:=floor(t/h)*h] %>% 
  .[,by=bin,media:=mean(value)] %>% 
  .[,by=bin,':='(inf=quantile(value,.05),sup=quantile(value,.95))] %>% 
  .[,!c('variable','value')] %>% 
  unique %>%
  .[,by=bin,esperanza:=(pnorm(bin+h)-pnorm(bin))/h] %>% 
  melt(id.vars=c('t','bin')) %>% 
  # .[,col:=fifelse(variable %in% c('inf','sup'),factor('inf',levels=c('media','inf','sup','esperanza')),variable )] %>% 
  .[,col:=fifelse(variable %in% c('inf','sup'),factor('inf',levels=levels(variable)),variable )] %>% 
  ggplot(aes(x=t,y=value,col=col,group=variable))+
  geom_line()
  

p1+
  stat_function(data=data.frame(c(-4,4)),fun=dnorm)+
  geom_line(data=cbind(t,estimaciones) %>% 
              melt(id.vars='t') %>% 
              .[,bin:=floor(t/h)*h] %>% 
              .[,by=bin,media:=mean(value)],
            aes(t,media,col='media'),show.legend=F)+
  geom_ribbon(data=cbind(t,estimaciones) %>% 
              melt(id.vars='t') %>% 
              .[,bin:=floor(t/h)*h] %>% 
              .[,by=bin,':='(inf=quantile(value,.05),sup=quantile(value,.95))],
            aes(x=t,ymin=inf,ymax=sup,col='media'),alpha=.5,show.legend=F)+
  geom_line(data=cbind(t,estimaciones) %>% 
              melt(id.vars='t') %>% 
              .[,bin:=floor(t/h)*h] %>% 
              .[,by=bin,value2:=mean(value)] %>% 
              .[,!c('variable','value')] %>% 
              unique %>% 
              .[,by=bin,esperanza:=(pnorm(bin+h)-pnorm(bin))/h],
            aes(t,esperanza,col='verdadera'),show.legend=F)


ggplot()+
  stat_function(data=data.frame(c(-4,4)),fun=dnorm)+
  geom_line(data=cbind(t,estimaciones) %>% 
              melt(id.vars='t') %>% 
              .[,bin:=floor(t/h)*h] %>% 
              .[,by=bin,media:=mean(value)],
            aes(t,media,col='media'),show.legend=F)+
  geom_ribbon(data=cbind(t,estimaciones) %>% 
                melt(id.vars='t') %>% 
                .[,bin:=floor(t/h)*h] %>% 
                .[,by=bin,':='(inf=quantile(value,.05),sup=quantile(value,.95))],
              aes(x=t,ymin=inf,ymax=sup,col='media'),alpha=.5,show.legend=F)+
  geom_line(data=cbind(t,estimaciones) %>% 
              melt(id.vars='t') %>% 
              .[,bin:=floor(t/h)*h] %>% 
              .[,by=bin,value2:=mean(value)] %>% 
              .[,!c('variable','value')] %>% 
              unique %>% 
              .[,by=bin,esperanza:=(pnorm(bin+h)-pnorm(bin))/h],
            aes(t,esperanza,col='verdadera'),show.legend=F)


# con otro h --------------------------------------------------------------

h <- 2
bins <- seq(xmin,xmax,h)

set.seed(pi)
estimaciones <- setDT(lapply(seq_len(750),function(foo) {
  x <- rnorm(n)
  cont <- table(cut(x,breaks = bins))/n
  f_est <- stepfun(bins[2:(length(bins)-1)],cont/h)
  f_est(t)
}))


(p1 <- cbind(t,estimaciones) %>% 
    melt(id.vars='t') %>% 
    ggplot()+
    geom_bin2d(aes(x=t,y=value,group=t,fill=..density..),show.legend=F,binwidth=c(.04,.02))+
    scale_fill_gradient(name='Densidad',trans=reverselog_trans(10))+
    stat_function(fun=dnorm,show.legend=F))

p1+
  stat_function(data=data.frame(c(-4,4)),fun=dnorm)+
  geom_line(data=cbind(t,estimaciones) %>% 
              melt(id.vars='t') %>% 
              .[,bin:=floor(t/h)*h] %>% 
              .[,by=bin,value2:=mean(value)],
            aes(t,value2,col='media'),show.legend=F)+
  geom_ribbon(data=cbind(t,estimaciones) %>% 
                melt(id.vars='t') %>% 
                .[,bin:=floor(t/h)*h] %>% 
                .[,by=bin,':='(inf=quantile(value,.05),sup=quantile(value,.95))],
              aes(x=t,ymin=inf,ymax=sup,col='media'),alpha=.5,show.legend=F)+
  geom_line(data=cbind(t,estimaciones) %>% 
              melt(id.vars='t') %>% 
              .[,bin:=floor(t/h)*h] %>% 
              .[,by=bin,value2:=mean(value)] %>% 
              .[,!c('variable','value')] %>% 
              unique %>% 
              .[,by=bin,esperanza:=(pnorm(bin+h)-pnorm(bin))/h],
            aes(t,esperanza,col='verdadera'),show.legend=F)

# ggplot()+
#   stat_function(data=data.frame(c(-4,4)),fun=dnorm)+
#   geom_line(data=cbind(t,estimaciones) %>% 
#               melt(id.vars='t') %>% 
#               .[,bin:=floor(t/h)*h] %>% 
#               .[,by=bin,value2:=mean(value)],aes(t,value2,col='media'),show.legend=F)+
#   geom_ribbon(data=cbind(t,estimaciones) %>% 
#                 melt(id.vars='t') %>% 
#                 .[,bin:=floor(t/h)*h] %>% 
#                 .[,by=bin,':='(inf=quantile(value,.05),sup=quantile(value,.95))],aes(x=t,ymin=inf,ymax=sup,fill='media'),alpha=.5,show.legend=F)


# h = 4 -------------------------------------------------------------------

h <- 4
bins <- seq(xmin,xmax,h)

set.seed(pi)
estimaciones <- setDT(lapply(seq_len(750),function(foo) {
  x <- rnorm(n)
  cont <- table(cut(x,breaks = bins))/n
  f_est <- stepfun(bins[2:(length(bins)-1)],cont/h)
  f_est(t)
}))


(p1 <- cbind(t,estimaciones) %>% 
    melt(id.vars='t') %>% 
    ggplot()+
    geom_bin2d(aes(x=t,y=value,group=t,fill=..density..),show.legend=F,binwidth=c(h,1/n/h))+
    scale_fill_gradient(name='Densidad',trans=reverselog_trans(10))+
    stat_function(fun=dnorm,show.legend = F))

p1+
  stat_function(data=data.frame(c(-4,4)),fun=dnorm)+
  geom_line(data=cbind(t,estimaciones) %>% 
              melt(id.vars='t') %>% 
              .[,bin:=floor(t/h)*h] %>% 
              .[,by=bin,value2:=mean(value)],
            aes(t,value2,col='media'),show.legend=F)+
  geom_ribbon(data=cbind(t,estimaciones) %>% 
                melt(id.vars='t') %>% 
                .[,bin:=floor(t/h)*h] %>% 
                .[,by=bin,':='(inf=quantile(value,.05),sup=quantile(value,.95))],
              aes(x=t,ymin=inf,ymax=sup,col='media'),alpha=.5,show.legend=F)+
  geom_line(data=cbind(t,estimaciones) %>% 
              melt(id.vars='t') %>% 
              .[,bin:=floor(t/h)*h] %>% 
              .[,by=bin,value2:=mean(value)] %>% 
              .[,!c('variable','value')] %>% 
              unique %>% 
              .[,by=bin,esperanza:=(pnorm(bin+h)-pnorm(bin))/h],
            aes(t,esperanza,col='verdadera'),show.legend=F)

# ggplot()+
#   stat_function(data=data.frame(c(-4,4)),fun=dnorm)+
#   geom_line(data=cbind(t,estimaciones) %>% 
#               melt(id.vars='t') %>% 
#               .[,bin:=floor(t/h)*h] %>% 
#               .[,by=bin,value2:=mean(value)],aes(t,value2,col='media'),show.legend=F)+
#   geom_ribbon(data=cbind(t,estimaciones) %>% 
#                 melt(id.vars='t') %>% 
#                 .[,bin:=floor(t/h)*h] %>% 
#                 .[,by=bin,':='(inf=quantile(value,.05),sup=quantile(value,.95))],aes(x=t,ymin=inf,ymax=sup,fill='media'),alpha=.5,show.legend=F)


# covarianzas entre las estimaciones en cada bin ???? --------------------------------------------------------


h <- .5
bins <- seq(xmin,xmax,h)

set.seed(pi)
estimaciones <- setDT(lapply(seq_len(1000),function(foo) {
  x <- rnorm(n)
  cont <- table(cut(x,breaks = bins))/n
  f_est <- stepfun(bins[2:(length(bins)-1)],cont/h)
  f_est(t)
}))

cbind(t,estimaciones,bin=floor(t/h)*h) %>%
  .[bin %in% unique(bin)[c(8,10)],by=bin,.SD[1]] %>%
  # .[bin %in% unique(bin)[c(80, 81)],by=bin,.SD[1]] %>% 
  .[,!c('bin')] %>%
  t %>% 
  as.data.table() %>% 
  {ggplot(data=.[-1],mapping = aes(V1,V2))+geom_point()+
      xlab(paste('Estimación en x =',.[1,1]))+ylab(paste('Estimación en x =',.[1,2]))}

# cbind(t,estimaciones) %>% 
#   melt(id.vars='t') %>% 
#   .[,bin:=floor(t/h)*h] %>% 
#   .[,.SD[1],by=.(bin,variable)] %>%
#   .[,value2:=ifelse(bin!=2,sum(.SD[bin!=2,value/.N]),value),by=variable] %>% 
#   .[]
