
# -------------------------------------------------------------------------
library(WaveletComp)
library(TSA)

set.seed(1234)


n <- 1e3
# t <- 1:n
# x <- sin(2*pi/50*t)
t <- 1:n*2*pi/100
tt <- c(300,500,750)
x <- c(sin(t[1:tt[1]]),
       sin(2*t[(tt[1]+1):tt[2]]),
       .5*sin(8*t[(tt[2]+1):tt[3]]) + .7*sin(2*t[(tt[2]+1):tt[3]]),
       sin(4*t[(tt[3]+1):n]))
y <- sin(4*t)
plot(t, x, type = "l")



my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 200,
                        dt = 1, dj = 1/100,
                        lowerPeriod = 6,
                        upperPeriod = 128,
                        make.pval = FALSE, n.sim = 10)

par(mfrow = c(2,1), mai = c(.3,.5,.3,.3))
plot(x, type = "l")
abline(v = t[aux[,1]], col = 2, lwd = 2)
wt.image(my.w, n.levels = 250, plot.legend = FALSE)


# -------------------------------------------------------------------------

library(ggplot2)
library(magrittr)
library(data.table)

par(mfrow = c(2,1), mai = c(.5,.5,.3,.3))
plot(t, x, type = "l", pch = 20)
abline(v = t[aux[,1]], col = 2, lwd = 2)
my.w$Ampl %>% t %>% image(z = ., x = t)
abline(v = t[aux[,1]], col = 2, lwd = 2)

spec.pgram(x, spans = c(2, 5))


plot(t, y, type = "l")
spec.pgram(y, spans = c(2, 5))

periodogram(x)
Mod(fft(x)^2) %>% plot(type = "l")

stats::spec(x)
stats::spec.ar(y)
spec(y)
tuneR::periodogram(y)


aux <- cbind(c(0, tt)+1, c(tt, n))
plot(t, x, type = "l")
abline(v = t[aux[,1]], col = 2, lwd = 4)
par(mfrow = c(4, 2), mai = c(.2,.2,.3,.3))
sapply(1:nrow(aux), \(j) {
  plot(t, x, type = "n")
  # abline(v = t[aux[,1]], col = 2, lwd = 2)
  rect(xleft = t[aux[j,1]], xright = t[aux[j,2]], ybottom = min(x), ytop = max(x), col = 2, density = 10)
  lines(t, x)
  spec(x[aux[j,1]:aux[j,2]])
})
