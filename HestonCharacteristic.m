function [y] = HestonCharacteristic(u,t,r,q,S0,hestonParameters)

kappa=hestonParameters(1);
eta=hestonParameters(2);
theta=hestonParameters(3);
corr=hestonParameters(4);
sig=hestonParameters(5);

d=((corr*theta*u*1i-kappa).^2-theta^2*(-1i.*u-u.^2)).^(1/2);
g=(kappa-corr*theta*u*1i-d)./(kappa-corr*theta*u*1i+d);
y=exp(1i.*u.*(log(S0)+(r-q).*t)).*exp(eta.*kappa.*theta^(-2).*((kappa-corr*theta*u*1i-d).*t-2*log((1-g.*exp(-d*t))./(1-g)))).*exp(sig^2*theta^(-2)*(kappa-corr*theta*u*1i-d).*(1-exp(-d.*t))./(1-g.*exp(-d.*t)));