test_fun <- function(t, y){
  dy <- vector(length = 2)
  dy[1] <- -sin(t)/((1+exp(2*t))^(1/2))+y[1]*(y[1]^2+y[2]^2-1)
  dy[2] <- cos(t)/((1+exp(2*t))^(1/2))+y[2]*(y[1]^2+y[2]^2-1)
  return (dy)
}

main_fun <- function(t, y){
  dy <- vector(length = 2)
  dy[1] <- g*(1-y[1])*y[1]-a*y[1]*y[2]/(b+y[1])
  dy[2] <- (1-y[2]/y[1])*y[2]
  return (dy)
}

RK3 <- function(f, tn, y0){
  y <- rep(y0, N+1)
  dim(y) <- c(2, N+1)
  for (i in 1:N){
    k1 <- f(tn[i], y[,i])
    k2 <- f(tn[i]+h/2, y[,i]+h/2*k1)
    k3 <- f(tn[i]+3/4*h, y[,i]+3/4*h*k2)
    y[,i+1] <- y[,i]+h*(2*k1+3*k2+4*k3)/9
  }
  return (y)
}

t1 <-  0
t2 <-  5
y0 = c(2^(1/2)/2, 0)
y <- RK3(test_fun, tn, y0)
df <- data.frame(x=cos(tn)/(1+exp(2*tn))^(1/2), y=sin(tn)/(1+exp(2*tn))^(1/2))
plot(df, type='l')
plot(x=y[1,], y=y[2,], type='l', xlab = 'y1', ylab = 'y2')

y0 = c(2^(1/2)/2, 0)
num_nodes <- (1:10)*10
hs <-  (t2-t1)/num_nodes
error <- c()
for (i in 1:length(hs)){
  h <- hs[i]
  N <- num_nodes[i]
  tn <-  t1 + (0:N)*h
  y <- RK3(test_fun, tn, y0)
  df <- data.frame(x=cos(tn)/(1+exp(2*tn))^(1/2), y=sin(tn)/(1+exp(2*tn))^(1/2))
  error <- c(error, max(c(abs(df$x-y[1,]), abs(df$y-y[2,]))))
}
plot(x = hs, y = error, type='l', xlab='Шаг', ylab = 'Погрешность')
plot(x = hs, y = error/hs^3, type='l', xlab='Шаг', ylab = 'Погрешность на шаг')

t1 <- 0
t2 <- 50
N <- 500
h <- (t2-t1)/N
tn <- t1 + (0:N)*h
y0 <- c(0.5, 0.5)
b <- 0.13
g <- 5
for (a in (1:30)){
  y <- RK3(main_fun, tn, y0)
  plot(x=y[1,], y = y[2,], type='l', main = y0, xlab = "X(t)", ylab = "Y(t)")
  plot(y[1,], type='l', main = a, xlab = "t", ylab = "X(t)")
  plot(y[2,], type='l', main = a, xlab = "t", ylab = "Y(t)")
}
b <- 0.13
a <- 10
for (g in (1:20)){
  y <- RK3(main_fun, tn, y0)
  plot(x=y[1,], y = y[2,], type='l', main = g, xlab = "X(t)", ylab = "Y(t)")
  plot(y[1,], type='l', main = a, xlab = "t", ylab = "X(t)")
  plot(y[2,], type='l', main = a, xlab = "t", ylab = "Y(t)")
}
g <- 12
a <- 10
for (b in (1:10)/10){
  y <- RK3(main_fun, tn, y0)
  plot(x=y[1,], y = y[2,], type='l', main = b, xlab = "X(t)", ylab = "Y(t)")
  plot(y[1,], type='l', main = a, xlab = "t", ylab = "X(t)")
  plot(y[2,], type='l', main = a, xlab = "t", ylab = "Y(t)")
}
