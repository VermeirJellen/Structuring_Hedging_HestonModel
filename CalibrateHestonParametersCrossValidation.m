function [hestonParameters outOfSampleRMSE outOfSampleRMSEAdjusted] = CalibrateHestonParametersCrossValidation(marketData,carrMadanPrecomputed,initialGuessHeston,weights,plotOption)

if(nargin < 2 | isempty(carrMadanPrecomputed)) % Compute Carr Madan Parameters, if missing
	N=4096; alpha=1.5; gridSpace=0.25; simpsonIntegrand=1;
	carrMadanPrecomputed = PrecomputationCarrMadanParameters(N,alpha,gridSpace,simpsonIntegrand); % Precompute modelindependent carr-madan parameters and variables
end;
if(nargin < 3 | isempty(initialGuessHeston)) % Declare initial Heston parameters, if missing
	kappa=7; eta=0.7; theta=0.7; corr=-0.5; sig=0.18; % Heston model, initial guess
	initialGuessHeston = [kappa eta theta corr sig];
end;
if(nargin < 4 | isempty(weights)) % compute weights, if missing
	weights = ones(1,length(marketData.T)); % use equal weighted approach
end;
if(nargin < 5)
	plotOption=false;
end;

% RANDOMLY SPLIT UP MARKETDATA IN TRAINING SET AND VALIDATION SET
percentageTest = 0.25; % Use 25% of the options for out of sample validation
trainingIndices = 1:length(marketData.T);
outOfSampleIndices = trainingIndices(randsample(trainingIndices,floor(length(trainingIndices)*percentageTest))); % 
trainingIndices(outOfSampleIndices) = []; % remove testsamples from training set

% TRAINING DATA
S0 = marketData.S0; q = marketData.q; T = marketData.T(trainingIndices); K = marketData.K(trainingIndices);
optionType = marketData.optionType(trainingIndices); optionMidMarket = marketData.optionMidMarket(trainingIndices);
bid = marketData.bid(trainingIndices); ask = marketData.ask(trainingIndices); r = marketData.r(trainingIndices);
weightsOptimization = weights(trainingIndices);

% TESTING DATA
TTesting = marketData.T(outOfSampleIndices); KTesting = marketData.K(outOfSampleIndices);
optionTypeTesting = marketData.optionType(outOfSampleIndices); optionMidMarketTesting = marketData.optionMidMarket(outOfSampleIndices);
bidTesting = marketData.bid(outOfSampleIndices); askTesting = marketData.ask(outOfSampleIndices); rTesting = marketData.r(outOfSampleIndices);


% PERFORM CROSS VALIDATION
tic;
nrIterationsBetweenValidation=75; converged=false; currentBestParameterGuess = initialGuessHeston; 
currentBestOutOfSampleRMSE = 10000; currentBestOutOfSampleRMSEAdjusted = 10000;
while(not(converged))

	newParameterGuess = HestonCalibrationLocalOptimization(currentBestParameterGuess,K,T,S0,r,q,optionType,carrMadanPrecomputed,optionMidMarket,weightsOptimization,false,nrIterationsBetweenValidation);
	
	newOutOfSamplePrices = HestonPricer(newParameterGuess,KTesting,TTesting,S0,rTesting,q,optionTypeTesting,carrMadanPrecomputed);
	newOutOfSampleRMSE = RMSE(newOutOfSamplePrices,optionMidMarketTesting);
	newOutOfSampleRMSEAdjusted = RMSESpreadAdjusted(newOutOfSamplePrices,optionMidMarketTesting,bidTesting,askTesting);
	
	if(newOutOfSampleRMSEAdjusted < currentBestOutOfSampleRMSEAdjusted)
		currentBestParameterGuess = newParameterGuess;
		currentBestOutOfSampleRMSE = newOutOfSampleRMSE;
		currentBestOutOfSampleRMSEAdjusted = newOutOfSampleRMSEAdjusted;
	else
		converged=true;
	end
end
toc;

hestonParameters = currentBestParameterGuess;
outOfSampleRMSE = currentBestOutOfSampleRMSE;
outOfSampleRMSEAdjusted  = currentBestOutOfSampleRMSEAdjusted;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%
% plot in sample and out of sample results
if(plotOption)
	inSamplePrices = HestonPricer(hestonParameters,K,T,S0,r,q,optionType,carrMadanPrecomputed);
	outOfSamplePrices = HestonPricer(hestonParameters,KTesting,TTesting,S0,rTesting,q,optionTypeTesting,carrMadanPrecomputed);
	
	figure;
	subplot(2,1,1);
	l=1:length(T); plot(l,optionMidMarket(l),'g.',l,bid(l),'r.',l,ask(l),'r.',l,inSamplePrices(l),'co');
	title('In sample results'); xlabel('option'); ylabel('price');
	legend('Mid Market','Bid','Ask','Heston Price','Location','northwest');
	
	subplot(2,1,2);
	l=1:length(TTesting); 
	plot(l,optionMidMarketTesting(l),'g.',l,bidTesting(l),'r.',l,askTesting(l),'r.',l,outOfSamplePrices(l),'co');
	title('Out Of Sample Results'); xlabel('option'); ylabel('price');
	legend('Mid Market','Bid','Ask','Heston Price','Location','northwest');
end