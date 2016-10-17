function minimizationCriterionBlackScholesExotic = MinimizationCriterionBlackScholesExotic(sigmaCalibration,S0,K,B,T,r,q,optionType,referencePrices)

pricesAndGreeks = BSBarrierPricerAndGreeksAnalytical(sigmaCalibration,S0,K,B,T,r,q,optionType); % Function currently only works for 1 option
minimizationCriterionBlackScholesExotic = sum((pricesAndGreeks(:,1)-referencePrices).^2);