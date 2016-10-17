% Pricer for one exotic DOBC option only
% For bulk matrix pricing, code can be re-implemented similarly to "MonteCarloSimulationRobustness"
% Last two parameters determine usage of antithetic variable and/or control variate
function [priceDOBC lowerLevel upperLevel] = MonteCarloPricerDOBCHeston(hestonParameters,K,T,S0,barrierPercentage,r,q,nrPeriods,nrPaths,confidenceLevel,useAntiThetic,useControlVariate)

B=barrierPercentage*S0;

if nargin < 11 | not(useAntiThetic) % Do not use antithetic random variables
	stock = SimulateStockPathsHeston(hestonParameters,S0,T,r,q,nrPeriods,nrPaths,[],[],false);
	payoffDOBC = max((min(stock,[],2)-B)./abs(min(stock,[],2)-B),0).*max((stock(:,nrPeriods+1)-K),0);
	nrSamples=nrPaths;
else % Use antithetic variable
	[stock, antiStock] = SimulateStockPathsHestonAntiThetic(hestonParameters,S0,T,r,q,nrPeriods,nrPaths,[],[],false);
	payoffDOBC=max((min(stock,[],2)-B)./abs(min(stock,[],2)-B),0).*max((stock(:,nrPeriods+1)-K),0);
	payoffDOBCNeg = max((min(antiStock,[],2)-B)./abs(min(antiStock,[],2)-B),0).*max((antiStock(:,nrPeriods+1)-K),0);
	nrSamples=nrPaths*2; % we have x2 the amount of individual samples
end;

if nargin==12 & useControlVariate % Use Stock price at maturity as control variate
 	N=nrPeriods; M=nrPaths;
	
	stMaturity = stock(:,N+1); stAverage = mean(stMaturity); averagePayoffDOBC = mean(payoffDOBC);
	A = (sum((payoffDOBC-averagePayoffDOBC).*(stMaturity-stAverage))); % Cov(DOBC,covariate) * nrSamples
	B = sum((stMaturity-stAverage).^2); % Var(Stock) * nrSamples
	alphaN = A/B; % alphaN = Cov(DOBC,stock) / Var(covariate)
	payoffDOBC = payoffDOBC - alphaN*(stMaturity-stAverage); % reduced variance average payoffs
	
	if useAntiThetic % Apply control variate logic to antithetic results
		stAntiMaturity = antiStock(:,N+1); stAntiAverage = mean(stAntiMaturity); averagePayoffDOBCNeg = mean(payoffDOBCNeg);
		ANeg = (sum((payoffDOBCNeg-averagePayoffDOBCNeg).*(stAntiMaturity-stAntiAverage)));
		BNeg = sum((stAntiMaturity-stAntiAverage).^2); 
		alphaNNeg = ANeg/BNeg;
		payoffDOBCNeg = payoffDOBCNeg - alphaNNeg*(stAntiMaturity-stAntiAverage);
	end
end

if(nargin >= 11 & useAntiThetic)
	 payoffDOBC = (payoffDOBC+payoffDOBCNeg)./2;
end;

% calculate final price
averagePayoffDOBC = mean(payoffDOBC);
priceDOBC = exp(-r*T)*averagePayoffDOBC;

% Calculate confidence levels
sampleVar = sum((payoffDOBC - averagePayoffDOBC).^2)/(length(payoffDOBC) - 1);
sampleStd = sqrt(sampleVar) / sqrt(nrSamples);
lowerLevel = exp(-r*T)*(averagePayoffDOBC + norminv(confidenceLevel/2)*sampleStd);
upperLevel = exp(-r*T)*(averagePayoffDOBC - norminv(confidenceLevel/2)*sampleStd);