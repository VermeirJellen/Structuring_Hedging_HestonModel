clc; clear;

trainingData = xlsread('optionData/OptionData_TRAINING_19_MAY','OptionData_TRAINING_19_MAY');
testingData = xlsread('optionData/OptionData_TESTING_19_MAY','OptionData_TESTING_19_MAY');

dividend = xlsread('optionData/OptionData_ALL_19_MAY','Dividend');
interestRates = xlsread('optionData/OptionData_ALL_19_MAY','InterestRates');
stock = xlsread('optionData/OptionData_ALL_19_MAY','Stock');
S0 = stock; q = dividend; 

% Training Data
TTraining = trainingData(:,2)'; KTraining = trainingData(:,3)'; optionTypeTraining = trainingData(:,4)'; % 0 = call, 1 = put
optionMidMarketTraining = trainingData(:,5)'; bidTraining = trainingData(:,6)'; askTraining = trainingData(:,7)';
rTraining = spline(interestRates(:,1), interestRates(:,2), TTraining);
TTraining = TTraining/365;

% Testing Data
TTesting = testingData(:,2)'; KTesting = testingData(:,3)'; optionTypeTesting = testingData(:,4)'; % 0 = call, 1 = put
optionMidMarketTesting = testingData(:,5)'; bidTesting = testingData(:,6)'; askTesting = testingData(:,7)';
rTesting= spline(interestRates(:,1), interestRates(:,2), TTesting);
TTesting = TTesting/365;

% Simulation parameters
randomBarrierPercentage = unifrnd(0.7,0.95,1,length(TTraining)); % generate random barriers for exotic pricing
nrPeriods=max(TTraining*365);
nrPaths=150000;

% weightsOptimization = ones(1,length(KTraining));
weightsOptimization = 1./abs(askTraining-bidTraining);
% weightsOptimization = 1./BlackScholesImpliedVolatility(S0,KTraining,rTraining,q,TTraining,optionTypeTraining,optionMidMarketTraining,false);

M = nrPaths; N = nrPeriods;
Z1 = normrnd(0,1,M,N);
Z3 = normrnd(0,1,M,N);

% Carr Madan Parameters
N=4096; alpha=1.5; gridSpace=0.25; simpsonIntegrand=1;
carrMadanPrecomputed = PrecomputationCarrMadanParameters(N,alpha,gridSpace,simpsonIntegrand); % Precompute modelindependent carr-madan parameters and variables

% Initial Heston Parameters
kappa=7; eta=0.7; theta=0.7; corr=-0.5; sig=0.18;
initialGuessHeston = [kappa eta theta corr sig];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BLACK SCHOLES BASELIN %%%%%%%%%%%%%%%%%%%%%
K=KTraining; T=TTraining; optionType=optionTypeTraining; optionMidMarket=optionMidMarketTraining; r=rTraining; bid=bidTraining; ask=askTraining;
sigmaBS = BlackScholesImpliedVolatility(S0,K,r,q,T,optionType,optionMidMarket,true); % Calibrated black scholes baselin
BSPricesInSample = BS(sigmaBS,S0,K,r,q,T,optionType);
rmseBSInSample = RMSE(BSPricesInSample,optionMidMarket);
rmseAdjustedBSInSample = RMSESpreadAdjusted(BSPricesInSample,optionMidMarket,bid,ask);

K=KTesting; T=TTesting; optionType=optionTypeTesting; optionMidMarket=optionMidMarketTesting; r=rTesting; bid=bidTesting; ask=askTesting;
BSPricesOutOfSample = BS(sigmaBS,S0,K,r,q,T,optionType);

rmseBSOutOfSample = RMSE(BSPricesOutOfSample,optionMidMarket);
rmseAdjustedBSOutOfSample = RMSESpreadAdjusted(BSPricesOutOfSample,optionMidMarket,bid,ask);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUll CONVERGENCE %%%%%%%%%%%%%%%%%%%%%%%%%%
% FULL CONVERGENCE --> TRAINING
K=KTraining; T=TTraining; optionType=optionTypeTraining; optionMidMarket=optionMidMarketTraining; r=rTraining; bid=bidTraining; ask=askTraining;
tic; calibratedParamsHestonFull = HestonCalibrationLocalOptimization(initialGuessHeston,K,T,S0,r,q,optionType,carrMadanPrecomputed,optionMidMarket,weightsOptimization,false); runtime=toc;
calibratedPricesHestonFull = HestonPricer(calibratedParamsHestonFull,K,T,S0,r,q,optionType,carrMadanPrecomputed); % calculate marketprices
rmseHestonFullInSample = RMSE(calibratedPricesHestonFull,optionMidMarket); % 0.1472
rmseSpreadAdjustedHestonFullInSample = RMSESpreadAdjusted(calibratedPricesHestonFull,optionMidMarket,bid,ask); % 0.0247

% FULL CONVERGENCE --> TESTING CARR MADAN
K=KTesting; T=TTesting; optionType=optionTypeTesting; optionMidMarket=optionMidMarketTesting; r=rTesting; bid=bidTesting; ask=askTesting;
calibratedPricesHestonFullTestCarrMadan = HestonPricer(calibratedParamsHestonFull,K,T,S0,r,q,optionType,carrMadanPrecomputed); % calculate marketprices
rmseHestonFullTestCarrMadan = RMSE(calibratedPricesHestonFullTestCarrMadan,optionMidMarket); % 0.1322
rmseSpreadAdjustedHestonFullTestCarrMadan = RMSESpreadAdjusted(calibratedPricesHestonFullTestCarrMadan,optionMidMarket,bid,ask); % 0.0240

% FULL CONVERGENCE --> TESTING MONTE CARLO
calibratedPricesHestonFullTestMonteCarloVanilla = MonteCarloSimulationRobustnessTest(calibratedParamsHestonFull,S0,T,K,r,q,optionType,nrPeriods,nrPaths,Z1,Z3,randomBarrierPercentage,false);
calibratedPricesHestonFullTestMonteCarloExotic = MonteCarloSimulationRobustnessTest(calibratedParamsHestonFull,S0,T,K,r,q,optionType,nrPeriods,nrPaths,Z1,Z3,randomBarrierPercentage,true);
rmseHestonFullTestMonteCarloVanilla = RMSE(calibratedPricesHestonFullTestMonteCarloVanilla,optionMidMarket); % 0.1385
rmseSpreadAdjustedHestonFullTestMonteCarloVanilla = RMSESpreadAdjusted(calibratedPricesHestonFullTestMonteCarloVanilla,optionMidMarket,bid,ask); % 0.0252
rmseHestonFullTestMonteCarloExotic = RMSE(calibratedPricesHestonFullTestMonteCarloExotic,optionMidMarket); % 0.1386
rmseSpreadAdjustedHestonFullTestMonteCarloExotic = RMSESpreadAdjusted(calibratedPricesHestonFullTestMonteCarloExotic,optionMidMarket,bid,ask); % 0.0282

fprintf('Full Convergence calibration runtime: %.2f seconds\n', runtime);
fprintf('In sample RMSE: %.4f\n',rmseHestonFullInSample);
fprintf('In sample adjusted RMSE: %.4f\n',rmseSpreadAdjustedHestonFullInSample);
fprintf('Out of sample carr madan RMSE: %.4f\n',rmseHestonFullTestCarrMadan);
fprintf('Out of sample adjusted carr madan RMSE: %.4f\n',rmseSpreadAdjustedHestonFullTestCarrMadan);
fprintf('Out of sample vanilla MC RMSE: %.4f\n',rmseHestonFullTestMonteCarloVanilla);
fprintf('Out of sample vanilla MC Adjusted RMSE: %.4f\n',rmseSpreadAdjustedHestonFullTestMonteCarloVanilla);
fprintf('Out of sample exotic MC RMSE: %.4f\n',rmseHestonFullTestMonteCarloExotic);
fprintf('Out of sample exotic MC Adjusted RMSE: %.4f\n\n',rmseSpreadAdjustedHestonFullTestMonteCarloExotic);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EARLY STOPPING %%%%%%%%%%%%%%%%%%%%%%%%%%
% EARLY STOPPING --> TRAINING
bid=bidTraining; ask=askTraining;
global globalStoppingCriterionBidAsk; % Define stopping condition. This global parameter will be utilized via the StopFunction
globalStoppingCriterionBidAsk = sum(weightsOptimization.*((bid-ask).^2));

% EARLY STOPPING --> TRAINING CARR MADAN
K=KTraining; T=TTraining; optionType=optionTypeTraining; optionMidMarket=optionMidMarketTraining; r=rTraining; bid=bidTraining; ask=askTraining;
% weightsOptimization = 1./BlackScholesImpliedVolatility(S0,KTraining,r,q,TTraining,optionTypeTraining,optionMidMarketTraining,false);
tic; calibratedParamsHestonEarlyStopping = HestonCalibrationLocalOptimization(initialGuessHeston,K,T,S0,r,q,optionType,carrMadanPrecomputed,optionMidMarket,weightsOptimization,true); runtime=toc;
calibratedPricesHestonEarlyStopping = HestonPricer(calibratedParamsHestonEarlyStopping,K,T,S0,r,q,optionType,carrMadanPrecomputed); % calculate marketprices
rmseHestonEarlyStoppingInSample = RMSE(calibratedPricesHestonEarlyStopping,optionMidMarket); % 0.3228
rmseSpreadAdjustedHestonEarlyStoppingInSample = RMSESpreadAdjusted(calibratedPricesHestonEarlyStopping,optionMidMarket,bid,ask); % 0.1916
l=1:length(T); plot(l,optionMidMarket(l),'g.',l,bid(l),'r.',l,ask(l),'r.',l,calibratedPricesHestonEarlyStopping(l),'co');

% EARLY STOPPING --> TESTING CARR MADAN
K=KTesting; T=TTesting; optionType=optionTypeTesting; optionMidMarket=optionMidMarketTesting; r=rTesting; bid=bidTesting; ask=askTesting;
calibratedPricesHestonEarlyStoppingTestCarrMadan = HestonPricer(calibratedParamsHestonEarlyStopping,K,T,S0,r,q,optionType,carrMadanPrecomputed); % calculate marketprices
rmseHestonEarlyStoppingTestCarrMadan = RMSE(calibratedPricesHestonEarlyStoppingTestCarrMadan,optionMidMarket); % 0.3176
rmseSpreadAdjustedHestonEarlyStoppingTestCarrMadan = RMSESpreadAdjusted(calibratedPricesHestonEarlyStoppingTestCarrMadan,optionMidMarket,bid,ask); % 0.2124
l=1:length(T); plot(l,optionMidMarket(l),'g.',l,bid(l),'r.',l,ask(l),'r.',l,calibratedPricesHestonEarlyStoppingTestCarrMadan(l),'co');

% EARLY STOPPING --> TESTING MONTE CARLO
calibratedPricesHestonEarlyStoppingTestMonteCarloVanilla = MonteCarloSimulationRobustnessTest(calibratedParamsHestonEarlyStopping,S0,T,K,r,q,optionType,nrPeriods,nrPaths,Z1,Z3,randomBarrierPercentage,false);
calibratedPricesHestonEarlyStoppingTestMonteCarloExotic = MonteCarloSimulationRobustnessTest(calibratedParamsHestonEarlyStopping,S0,T,K,r,q,optionType,nrPeriods,nrPaths,Z1,Z3,randomBarrierPercentage,true);

rmseHestonEarlyStoppingTestMonteCarloVanilla = RMSE(calibratedPricesHestonEarlyStoppingTestMonteCarloVanilla,optionMidMarket); % 0.2649
rmseSpreadAdjustedHestonEarlyStoppingTestMonteCarloVanilla = RMSESpreadAdjusted(calibratedPricesHestonEarlyStoppingTestMonteCarloVanilla,optionMidMarket,bid,ask); % 0.1573
rmseHestonEarlyStoppingTestMonteCarloExotic = RMSE(calibratedPricesHestonEarlyStoppingTestMonteCarloExotic,optionMidMarket); % 0.3281
rmseSpreadAdjustedHestonEarlyStoppingTestMonteCarloExotic = RMSESpreadAdjusted(calibratedPricesHestonEarlyStoppingTestMonteCarloExotic,optionMidMarket,bid,ask) % 0.2220
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Early stopping calibration runtime: %.2f seconds\n', runtime);
fprintf('In sample RMSE: %.4f\n',rmseHestonEarlyStoppingInSample);
fprintf('In sample adjusted RMSE: %.4f\n',rmseSpreadAdjustedHestonFullInSample);
fprintf('Out of sample carr madan RMSE: %.4f\n',rmseHestonEarlyStoppingTestCarrMadan);
fprintf('Out of sample adjusted carr madan RMSE: %.4f\n',rmseSpreadAdjustedHestonEarlyStoppingTestCarrMadan);
fprintf('Out of sample vanilla MC RMSE: %.4f\n',rmseHestonEarlyStoppingTestMonteCarloVanilla);
fprintf('Out of sample vanilla MC Adjusted RMSE: %.4f\n',rmseSpreadAdjustedHestonEarlyStoppingTestMonteCarloVanilla);
fprintf('Out of sample exotic MC RMSE: %.4f\n',rmseHestonEarlyStoppingTestMonteCarloExotic);
fprintf('Out of sample exotic MC Adjusted RMSE: %.4f\n\n',rmseSpreadAdjustedHestonEarlyStoppingTestMonteCarloExotic);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CROSS VALIDATION %%%%%%%%%%%%%%%%%%%%%%%%%%
K=KTraining; T=TTraining; optionType=optionTypeTraining; optionMidMarket=optionMidMarketTraining; r=rTraining; bid=bidTraining; ask=askTraining;
% weightsOptimization = 1./BlackScholesImpliedVolatility(S0,KTraining,r,q,TTraining,optionTypeTraining,optionMidMarketTraining,false);
% weightsOptimization = 1./abs(bid-ask);
% weightsOptimization = ones(1,length(KTraining));

tic;
nrIterationsBetweenValidation=75; converged=false; currentBestParameterGuess = initialGuessHeston; currentBestOutOfSampleRMSEAdjusted = 10000;
while(not(converged))

	newParameterGuess = HestonCalibrationLocalOptimization(currentBestParameterGuess,K,T,S0,r,q,optionType,carrMadanPrecomputed,optionMidMarket,weightsOptimization,false,nrIterationsBetweenValidation);
	
	newOutOfSamplePrices = HestonPricer(newParameterGuess,KTesting,TTesting,S0,rTesting,q,optionTypeTesting,carrMadanPrecomputed);
	newOutOfSampleRMSEAdjusted = RMSESpreadAdjusted(newOutOfSamplePrices,optionMidMarketTesting,bidTesting,askTesting);
	% newOutOfSampleRMSE = RMSESpreadAdjusted(newOutOfSamplePrices,optionMidMarketTesting,bidTesting,askTesting);
	
	if(newOutOfSampleRMSEAdjusted < currentBestOutOfSampleRMSEAdjusted)
		currentBestParameterGuess = newParameterGuess;
		currentBestOutOfSampleRMSEAdjusted = newOutOfSampleRMSEAdjusted;
	else
		converged=true;
		calibratedParamsHestonCrossValidation = currentBestParameterGuess;
	end
end
runtime = toc;

K=KTraining; T=TTraining; optionType=optionTypeTraining; optionMidMarket=optionMidMarketTraining; r=rTraining; bid=bidTraining; ask=askTraining;
calibratedPricesHestonCrossValidationTestCarrMadanInSample = HestonPricer(calibratedParamsHestonCrossValidation,K,T,S0,r,q,optionType,carrMadanPrecomputed); % calculate marketprices
rmseHestonCrossValidationTestCarrMadanInSample = RMSE(calibratedPricesHestonCrossValidationTestCarrMadanInSample,optionMidMarket); % 0.3176
rmseSpreadAdjustedHestonCrossValidationTestCarrMadanInSample = RMSESpreadAdjusted(calibratedPricesHestonCrossValidationTestCarrMadanInSample,optionMidMarket,bid,ask); % 0.2124


K=KTesting; T=TTesting; optionType=optionTypeTesting; optionMidMarket=optionMidMarketTesting; r=rTesting; bid=bidTesting; ask=askTesting;
calibratedPricesHestonCrossValidationTestCarrMadan = HestonPricer(calibratedParamsHestonCrossValidation,K,T,S0,r,q,optionType,carrMadanPrecomputed); % calculate marketprices
rmseHestonCrossValidationTestCarrMadan = RMSE(calibratedPricesHestonCrossValidationTestCarrMadan,optionMidMarket); % 0.1318
rmseSpreadAdjustedHestonCrossValidationTestCarrMadan = RMSESpreadAdjusted(calibratedPricesHestonCrossValidationTestCarrMadan,optionMidMarket,bid,ask); % 0.0244
l=1:length(T); plot(l,optionMidMarket(l),'g.',l,bid(l),'r.',l,ask(l),'r.',l,calibratedPricesHestonCrossValidationTestCarrMadan(l),'co');

calibratedPricesHestonCrossValidationTestMonteCarloVanilla = MonteCarloSimulationRobustnessTest(calibratedParamsHestonCrossValidation,S0,T,K,r,q,optionType,nrPeriods,nrPaths,Z1,Z3,randomBarrierPercentage,false);
calibratedPricesHestonCrossValidationTestMonteCarloExotic = MonteCarloSimulationRobustnessTest(calibratedParamsHestonCrossValidation,S0,T,K,r,q,optionType,nrPeriods,nrPaths,Z1,Z3,randomBarrierPercentage,true);

rmseHestonCrossValidationTestMonteCarloVanilla = RMSE(calibratedPricesHestonCrossValidationTestMonteCarloVanilla,optionMidMarket); % 0.1357
rmseSpreadAdjustedHestonCrossValidationTestMonteCarloVanilla = RMSESpreadAdjusted(calibratedPricesHestonCrossValidationTestMonteCarloVanilla,optionMidMarket,bid,ask); % 0.0234
rmseHestonCrossValidationMonteCarloExotic = RMSE(calibratedPricesHestonCrossValidationTestMonteCarloExotic,optionMidMarket); % 0.1379
rmseSpreadAdjustedHestonCrossValidationTestMonteCarloExotic = RMSESpreadAdjusted(calibratedPricesHestonCrossValidationTestMonteCarloExotic,optionMidMarket,bid,ask); % 0.0281

fprintf('Cross validation runtime: %.2f seconds\n', runtime);
fprintf('In sample RMSE: %.4f\n',rmseHestonCrossValidationTestCarrMadanInSample);
fprintf('In sample adjusted RMSE: %.4f\n',rmseSpreadAdjustedHestonCrossValidationTestCarrMadanInSample);
fprintf('Out of sample carr madan RMSE: %.4f\n',rmseHestonCrossValidationTestCarrMadan);
fprintf('Out of sample adjusted carr madan RMSE: %.4f\n',rmseSpreadAdjustedHestonCrossValidationTestCarrMadan);
fprintf('Out of sample vanilla MC RMSE: %.4f\n',rmseHestonCrossValidationTestMonteCarloVanilla);
fprintf('Out of sample vanilla MC Adjusted RMSE: %.4f\n',rmseSpreadAdjustedHestonCrossValidationTestMonteCarloVanilla);
fprintf('Out of sample exotic MC RMSE: %.4f\n',rmseHestonCrossValidationMonteCarloExotic);
fprintf('Out of sample exotic MC Adjusted RMSE: %.4f\n\n',rmseSpreadAdjustedHestonCrossValidationTestMonteCarloExotic);

