% Perform carr madan for the options that mature at timestamp T
function [y] = CarrMadanHestonOptimized(K,T,S0,r,q,type,cmPrecomputed,charPrecomputed)

rho = exp(-r.*T).*HestonCharacteristicOptimized(cmPrecomputed.u,T,r,q,charPrecomputed)./cmPrecomputed.rhoDenominator;  
a = real(fft(rho.*cmPrecomputed.fftMultiplier, cmPrecomputed.N));
CallPrices = (1/pi)*exp(-cmPrecomputed.alpha*cmPrecomputed.k).*a;       

y = spline(cmPrecomputed.realStrikes,CallPrices,K);
y(type==1) = y(type==1) + K(type==1)*exp(-r*T)- exp(-q*T)*S0; 