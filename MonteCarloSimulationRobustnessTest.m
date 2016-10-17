% Function can be utilized to verify out of sample robustness of the calibration for pricing exotics
function prices = MonteCarloSimulationRobustnessTest(calibratedParams,S0,T,K,r,q,optionType,nrPeriods,nrPaths,Z1,Z3,randomBarriers,useExotics)

results = zeros(1,length(T));
timeStamps = unique(T);

for(i=1:length(timeStamps))
	logicalVector = (T==timeStamps(i));
	KK = K(logicalVector);
	B=randomBarriers(logicalVector)*S0;
	optionTypeCm = optionType(logicalVector);
	
	first=find(logicalVector,1);
	tCm=T(first);
	rCm=r(first);

	if(useExotics) % use DOBC and DIBC to price call option
		Stock = SimulateStockPathsHeston(calibratedParams,S0,tCm,rCm,q,nrPeriods,nrPaths,Z1,Z3,false); % Simulate stock prices (return full paths, no antithetic var, no control variate)	
		payoffMatrixDOBC = zeros(length(Stock),length(KK));
		payoffMatrixDIBC = payoffMatrixDOBC;
		for (j=1:length(KK)) % 1 column contains payoff of each option KK
			payoffMatrixDOBC(:,j) = max((min(Stock,[],2)-B(j))./abs(min(Stock,[],2)-B(j)),0).*max((Stock(:,nrPeriods+1)-KK(j))*exp(-rCm*tCm),0);
			payoffMatrixDIBC(:,j) = max((B(j)-min(Stock,[],2))./abs(B(j)-min(Stock,[],2)),0).*max((Stock(:,nrPeriods+1)-KK(j))*exp(-rCm*tCm),0);
		end;
		mcPrice = mean(payoffMatrixDOBC) + mean(payoffMatrixDIBC); % price of call = price DOBC + price DIBC
	else % simulate call prices directly
		Stock = SimulateStockPathsHeston(calibratedParams,S0,tCm,rCm,q,nrPeriods,nrPaths,Z1,Z3,true); % Simulate stock prices (return full paths, no antithetic var, no control variate)	
		payoffMatrix = zeros(length(Stock),length(KK));
		for (j=1:length(KK)) % 1 column contains discounted payoff of each option KK
			payoffMatrix(:,j) = max(Stock-KK(j)*exp(-rCm*tCm),0);
		end;
		mcPrice = mean(payoffMatrix);
	end;
	
	mcPrice(optionTypeCm==1) = mcPrice(optionTypeCm==1) + KK(optionTypeCm==1)*exp(-rCm.*tCm)-S0*exp(-q*tCm); % Put Call parity
	results(logicalVector) = mcPrice;
end;

prices=results;