function out = Black_scholes(sigma, S0, K, r, q, T) 

d_1 = (log(S0./K)+(r - q +sigma^2/2).*T)./(sigma*sqrt(T));
d_2 = d_1 - sigma.*sqrt(T);
out = exp(-q.*T)*S0.*normcdf(d_1)-K.*exp(-r.*T).*normcdf(d_2);
