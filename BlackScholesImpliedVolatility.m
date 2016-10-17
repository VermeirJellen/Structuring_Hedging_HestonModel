function [impliedVolatility] = BlackScholesImpliedVolatility(S0,K,r,q,T,optionType,price,calibrateAllData)

optimizationSettings=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals',10000);

if(calibrateAllData) % find implied volatility over the whole dataset
	impliedVolatility = fminsearch(@(sigma)MinimizationCriterionBlackScholes(sigma,S0,K,r,q,T,optionType,price),0.5,optimizationSettings);
else  %% find the implied volatility for each option
	impliedVolatility = zeros(1,length(T));
	for i=1:length(T)
		impliedVolatility(1,i)= fminsearch(@(sigma)MinimizationCriterionBlackScholes(sigma,S0,K(i),r(i),q,T(i),optionType(i),price(i)),0.5,optimizationSettings);
		% impliedVolatilityBuiltIn(1,i) = blsimpv(S0,K(i),r(i),T(i),abs(price(i)),[],q,[],not(optionType(i)));
	end;
end;