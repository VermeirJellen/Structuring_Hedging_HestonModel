function out = HestonPricerWithoutPrecomputation(hestonParameters,K,T,S0,r,q,optionType)
		
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

	results(logicalVector) = CarrMadanHeston(KK,tCm,S0,rCm,q,optionTypeCm,hestonParameters);
end

out=results;