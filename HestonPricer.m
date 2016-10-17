function out = HestonPricer(hestonParameters,K,T,S0,r,q,optionType,carrMadanPrecomputed)

if nargin < 8
	N=4096; alpha=1.5; gridSpace=0.25; simpsonIntegrand=1;
	carrMadanPrecomputed = PrecomputationCarrMadanParameters(N,alpha,gridSpace,simpsonIntegrand);
end			
characteristicFunctionPrecomputed = PrecomputationHestonCharacteristicFunctionParameters(hestonParameters,carrMadanPrecomputed.u,S0);
			
timeStamps = unique(T);
results = zeros(1,length(T));

% Iterate over the unique timeStamps. Process all options for identical timestamps simultaneously
for(i=1:length(timeStamps))
	logicalVector = (T==timeStamps(i));
	KK = K(logicalVector);
	optionTypeCm = optionType(logicalVector);
	
	first=find(logicalVector,1);
	tCm=T(first);
	rCm=r(first);

	results(logicalVector) = CarrMadanHestonOptimized(KK,tCm,S0,rCm,q,optionTypeCm,carrMadanPrecomputed,characteristicFunctionPrecomputed);
end

out=results;