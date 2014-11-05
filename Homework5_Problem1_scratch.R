library(ggplot2)

LBandFreq  = 1.5e9 # Hz
SBandFreq  = 3e9 # Hz
XBandFreq  = 10e9 # Hz
KuBandFreq = 15e9 # Hz

c = 3e8 # m/s

lambda = c / LBandFreq
#lambda = 0.1

i = 2
N = 4*i
b0 = pi/2
k = (2*pi)/lambda
d = lambda/2


psi = function(theta, alpha, k, d) {
  if(length(alpha) == 1) alpha = rep(alpha, length(theta))
  
  if(length(theta) != length(alpha)) {
    cat('The theta and the alpha arrays are not of the same length\n')
    return(NA)
  }
  
  ret = k*d*cos(theta) + alpha
  return(ret)
}

f = function(theta, alpha, k, d, n, c) {
  if(length(alpha) == 1) alpha = rep(alpha, n)
  if(length(alpha) != n) {
    cat('The stearing phase is not the same length as the number of elements\n')
    return(NA)
  }
  
  p = sapply(1:n, FUN = function(x) psi(theta, alpha[x], k, d)) 
#  p = psi(theta, alpha, k, d)
  p = matrix(p, ncol = length(theta), byrow = TRUE)
  i = complex(imaginary = 1)
  ns = (1:n)-1
  ns = matrix(rep(ns, length(theta)), ncol = length(theta), byrow = FALSE)
  cs = matrix(rep(c, length(theta)), ncol = length(theta), byrow = FALSE)
  terms = cs*exp(ns*i*p)
  ret = colSums(terms)  
  ret = Re(sqrt(ret * Conj(ret)))
  return(ret)
}

f_ana = function(theta, alpha, k, d, n, c) {
  p = psi(theta, alpha, k, d)
  A = sin(n*p/2)
  B = sin((n*p/2)*(5/8))
  C = sin((n*p/2)*(3/8))
  ret = sqrt((1/(sin(p/2)*sin(p/2)))*((A*A)+(B*B)+(C*C) - (2*A*B*cos((n*p/2)*(3/8))) + (2*A*C*cos((n*p/2)*(5/8))) - (2*B*C*cos((n*p/2)*(2/8)))))
  ret = ret*(1/((n/8)*sqrt(98)))
  return(ret)
}

f_ana2 = function(theta, alpha, k, d, n, c) {
  i = complex(imaginary = 1)
  p = psi(theta, alpha, k, d)
  A = sin(n*p/2)
  B = sin((n*p/2)*(5/8))
  C = sin((n*p/2)*(3/8))
  D = sin(p/2)
  
  ret = (exp(i*(n-1)*p/2)*(A/D)) - (exp(i*((5*n/8)-1)*p/2)*(B/D)) + (exp(i*((3*n/8)-1)*p/2)*(C/D))
#  ret = (exp(i*(n-1)*p/2)*(A/D))
  ret = (1/((n/8)*sqrt(98)))*Re(sqrt(ret*Conj(ret)))
  
  return(ret)
}

C = rep(1, N)
# Remove center quarter elements
lo = (3*N)/8
hi = (5*N)/8
C[c(4,5)] = 0

alpha = k*d*(1 - ((2*b0)/pi))
alpha = rep(alpha, N)
#alpha[seq(1,N, by = 2)] = 2

theta = seq(0, pi, length.out = 1000)
res = f(theta, alpha, k, d, N, C)
#res = f_ana(theta, alpha[1], k, d, N, C)
resNorm = res / max(res)
res2 = f_ana(theta, alpha[1], k, d, N, C)
res2Norm = res2
#res2Norm = res2 / max(res2)

pdata = data.frame(theta = theta, cosTheta = cos(theta), sinTheta = sin(theta), psi = psi(theta, alpha, k, d), f = res, fN = resNorm, f2 = res2, fN2 = res2Norm)
rp = ggplot(pdata) + theme_bw() +
  geom_line(aes(x = cosTheta, y = fN), size = 1.5) +
  geom_line(aes(x = cosTheta, y = fN2), color = 'red')
rp




