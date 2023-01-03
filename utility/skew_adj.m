function sk = skew_adj(t,x)

% sk = skew_adj(t,x)
% Given the vector x (signal) and the scalar t, it computes the skew of x./(1 -t*x)
% JPM, April 27 2021, IO-CSIC, Aranjuez

xt = x./(1-t*x);
mu1 = mean(xt);
mu2 = mean((xt-mu1).^2);
sk = mean(((xt-mu1)/sqrt(mu2)).^3);