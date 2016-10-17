% Todo: Modify in such a way that sigma calibration can occur in mass
% Currently only works for non-vectorized inputs, BSBarrierPricerAndGreeksAnalytical should be modified to change this.
function [impliedVolatility] = BlackScholesImpliedVolatilityExotic(S0,K,B,T,r,q,optionType,price)

optimizationSettings=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals',10000);
impliedVolatility = fminsearch(@(sigma)MinimizationCriterionBlackScholesExotic(sigma,S0,K,B,T,r,q,optionType,price),0.5,optimizationSettings);

% if(calibrateAllData) % find implied volatility over the whole dataset
% 	impliedVolatility = fminsearch(@(sigma)MinimizationCriterionBlackScholesExotic(sigma,S0,K,B,T,r,q,optionType,price),0.5,optimizationSettings);
% else  %% find the implied volatility for each option
% 	impliedVolatility = zeros(1,length(T));
%	for i=1:length(T)
%		impliedVolatility(1,i)= fminsearch(@(sigma)MinimizationCriterionBlackScholesExotic(sigma,S0,K(i),B(i),T(i),r(i),q,optionType(i),price(i)),0.5,optimizationSettings);
%	end;
% end;