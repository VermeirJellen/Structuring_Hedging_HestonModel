clc; clear;

% CSV data generated via R script. 
option = xlsread('optionData/OptionData_ALL_19_MAY','OptionData_ALL_19_MAY');
dividend = xlsread('optionData/OptionData_ALL_19_MAY','Dividend');
interestRates = xlsread('optionData/OptionData_ALL_19_MAY','InterestRates');
stock = xlsread('optionData/OptionData_ALL_19_MAY','Stock');

S0 = stock; q = dividend; T = option(:,2)'; K = option(:,3)';
optionType = option(:,4)'; % 0 = call, 1 = put
optionMidMarket = option(:,5)'; bid = option(:,6)'; ask = option(:,7)';
r = spline(interestRates(:,1), interestRates(:,2), T);
T = T/365; % maturity in years
					
marketData = struct('S0',S0,'q',q,'T',T,'K',K,'optionType',optionType,'optionMidMarket',optionMidMarket,'bid',bid,'ask',ask,'r',r);


%% BLACK SCHOLES MODEL CALIBRATION %% 
% We will first calibrate the black Scholes model. This allows us to make a comparison between the black scholes model and the heston model
calibratedSigmaBlackScholes = BlackScholesImpliedVolatility(S0,K,r,q,T,optionType,optionMidMarket,true); % Note: true flag at the end indicates sigma calibration for whole data set
calibratedPricesBlackScholes = BS(calibratedSigmaBlackScholes,S0,K,r,q,T,optionType);
rmseBS = RMSE(calibratedPricesBlackScholes,optionMidMarket); % 0.2942
rmseAdjustedBS = RMSESpreadAdjusted(calibratedPricesBlackScholes,optionMidMarket,bid,ask); % 0.1717
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% HESTON MODEL CALIBRATION (CARR MADAN WITH CROSS VALIDATION APPROACH %%
% Carr Madan Parameters
N=4096; alpha=1.5; gridSpace=0.25; simpsonIntegrand=1;
carrMadanPrecomputed = PrecomputationCarrMadanParameters(N,alpha,gridSpace,simpsonIntegrand); % Precompute modelindependent carr-madan parameters and variables

% Initial educated guess for HestonParameters
kappa=7; eta=0.7; theta=0.7; corr=-0.5; sig=0.18;
initialGuessHeston = [kappa eta theta corr sig];

% 1. Equal weighting calibration (default)
% [hestonParametersEqual outOfSampleRMSE1 outOfSampleRMSEAdjusted1] = CalibrateHestonParametersCrossValidation(marketData,carrMadanPrecomputed,initialGuessHeston); % (17.1 seconds)
% outOfSampleRMSE1 = 0.1706 % outOfSampleAdjustedRMSE1 = 0.0248

% 2. Use bid ask spread as weights
% weightsOptimization = 1./abs(bid-ask);
% [hestonParametersSpread outOfSampleRMSE2 outOfSampleRMSEAdjusted2] = CalibrateHestonParametersCrossValidation(marketData,carrMadanPrecomputed,initialGuessHeston,weightsOptimization); % (20.3 seconds)
% outOfSampleRMSE2 = 0.1735 % outOfSampleAdjustedRMSE2 = 0.0332

% 3. Use Implied volatility as weights
weightsOptimization = 1./BlackScholesImpliedVolatility(S0,K,r,q,T,optionType,optionMidMarket,false);
[hestonParameters, outOfSampleRMSEHeston, outOfSampleRMSEAdjustedHeston] = CalibrateHestonParametersCrossValidation(marketData,carrMadanPrecomputed,initialGuessHeston,weightsOptimization,true); % (14.3 seconds)
% outOfSampleRMSE1 = 0.1367 % outOfSampleAdjustedRMSE1 = 0.0163
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outOfSampleRMSEHeston;
outOfSampleRMSEAdjustedHeston;

%% MONTE CARLO SIMULATION (No variance reduction, Anti thetic variable, Control variate, Anti thetic + control variate) 
B=0.70; KMaturity=S0; nrDays = 652; % 652 days until maturity
rMaturity=spline(interestRates(:,1), interestRates(:,2), nrDays); TMaturity=nrDays/365;
confidence=0.05; % parameter for confidence intervals

% No antithetic variable, no controlvariate, 150000 individual samples
nrPeriods=nrDays; nrPaths=150000;
tic; [priceDOBCNoVarReduction, lowerLevel, upperLevel] = MonteCarloPricerDOBCHeston(hestonParameters,KMaturity,TMaturity,S0,B,rMaturity,q,nrPeriods,nrPaths,confidence); runtime = toc; % 6.3 seconds
confidenceInterval = upperLevel-lowerLevel; % confidence interval (upperLevel-lowerLevel) =  0.0640
fprintf('MC without variance reduction price = %.2fUSD\n',priceDOBCNoVarReduction);
fprintf('MC without variance reduction confidence interval for %d paths and %d periods is %.4f\n',nrPaths,nrPeriods,confidenceInterval);
fprintf('MC without variance reduction runtime = %.2f seconds\n',runtime);

% Antithetic variable, no controlVariate, 75000 pairs
nrPeriods=nrDays; nrPaths=75000; % use only half the number of paths!
tic; [priceDOBCAntiThetic, lowerLevel, upperLevel] = MonteCarloPricerDOBCHeston(hestonParameters,KMaturity,TMaturity,S0,B,rMaturity,q,nrPeriods,nrPaths,confidence,true); runtime = toc; % 3.6 seconds
confidenceInterval = upperLevel-lowerLevel; % confidence interval (upperLevel-lowerLevel) =  0.0331
fprintf('MC with anti thetic random variables price = %.2fUSD\n',priceDOBCAntiThetic);
fprintf('MC with anti thetic random variables confidence interval for %d paths and %d periods is %.4f\n',nrPaths,nrPeriods,confidenceInterval);
fprintf('MC with anti thetic random variables runtime = %.2f seconds\n',runtime);

% Control variate, no antithetic variable, 150000 individual samples
nrPeriods=nrDays; nrPaths=150000; % use only half the number of paths!
tic; [priceDOBCControlVariate, lowerLevel, upperLevel] = MonteCarloPricerDOBCHeston(hestonParameters,KMaturity,TMaturity,S0,B,rMaturity,q,nrPeriods,nrPaths,confidence,false,true); runtime = toc; % 6.354003 seconds
confidenceInterval = upperLevel-lowerLevel; % confidence interval (upperLevel-lowerLevel) =  0.0367
fprintf('MC with control variate price = %.2fUSD\n',priceDOBCControlVariate);
fprintf('MC with control variate confidence interval for %d paths and %d periods is %.4f\n',nrPaths,nrPeriods,confidenceInterval);
fprintf('MC with control variate runtime = %.2f seconds\n',runtime);

% Control variate AND anti thetic variable, 75000 pairs
nrPeriods=nrDays; nrPaths=75000; % use only half the number of paths!
tic; [priceDOBCAntiTheticControlVariate, lowerLevel, upperLevel] = MonteCarloPricerDOBCHeston(hestonParameters,KMaturity,TMaturity,S0,B,rMaturity,q,nrPeriods,nrPaths,confidence,true,true); runtime = toc; % 3.6 seconds
confidenceInterval = upperLevel-lowerLevel; % confidence interval (upperLevel-lowerLevel) =  0.0367
fprintf('MC with control variate and anti thetic random variables price = %.2fUSD\n',priceDOBCAntiTheticControlVariate);
fprintf('MC with control variate and anti thetic random variables confidence interval for %d paths and %d periods is %.4f\n',nrPaths,nrPeriods,confidenceInterval);
fprintf('MC with control variate and anti thetic random variables runtime = %.2f seconds\n',runtime);



%% Intuition and some numbers that give insight on the relatively high DOBC price and vice versa low DIBC price %%
hestonPriceDOBC = priceDOBCAntiTheticControlVariate;
barrierLevel=0.70;
impliedVol = BlackScholesImpliedVolatilityExotic(S0,KMaturity,barrierLevel*S0,TMaturity,rMaturity,q,2,hestonPriceDOBC); % Implied volatility for DOBC matching heston price
priceAndGreeksDOBC = BSBarrierPricerAndGreeksAnalytical(impliedVol,S0,KMaturity,barrierLevel*S0,TMaturity,rMaturity,q,2);
priceAndGreeksDIBC = BSBarrierPricerAndGreeksAnalytical(impliedVol,S0,KMaturity,barrierLevel*S0,TMaturity,rMaturity,q,3);
priceAndGreeksCall = BSBarrierPricerAndGreeksAnalytical(impliedVol,S0,KMaturity,barrierLevel*S0,TMaturity,rMaturity,q,0);

fprintf('\nprice of the DOBC that we are interested in is %.2f\n',priceAndGreeksDOBC(1));
fprintf('In comparison, price of the corresponding DIBC with same T, B and K is only %.2f\n\n',priceAndGreeksDIBC(1));
fprintf('We can explain the low DIBC price by performing some stock simulations\n');

barrierLevel = 0.7; nrSimulations = 100000;
stockPaths = SimulateStockPathsHeston(hestonParameters,S0,TMaturity,rMaturity,q,nrDays,nrSimulations,[],[],false);
barrierBreachedMaturityLevels = stockPaths(find(min(stockPaths,[],2) < barrierLevel*S0), size(stockPaths,2)); % Stock price at maturity for paths that went below the barrier during the lifetime
barrierBreachedAndRecoveredMaturityLevels = barrierBreachedMaturityLevels(find(barrierBreachedMaturityLevels > S0)); % % Stock price at maturity for paths that went below the barrier, but finished in the money

nrBarrierBreached = length(barrierBreachedMaturityLevels); % There are 26391 out of 100000 stocks for which the barrier was breached during [0, T]
nrBarrierBreachedAndRecovered = length(barrierBreachedAndRecoveredMaturityLevels); % There were only 767 out of these 26391 stockpaths that recovered and finished in the money (>S0) at maturity
averagePercentageReturnBreachedAndRecovered = (mean(barrierBreachedAndRecoveredMaturityLevels)-S0)/S0; % The average (stock) profit of these 767 paths was only 7%

fprintf('Result of %d stocksimulations:\n',nrSimulations);
fprintf('%d breached the barrierLevel of %d percent\n',nrBarrierBreached,barrierLevel*100);
fprintf('%d of these "barrier-breached simulations" recovered and ended above the initial index price (St > S0)\n',nrBarrierBreachedAndRecovered);
fprintf('The final profit at maturity of these specific stock paths was on average %.2f percent\n',averagePercentageReturnBreachedAndRecovered*100);

fprintf('Only these specific stoch paths give a positive payoff for the DIBC\n');
fprintf('We can see that on average of all the stockpaths the expected payoff will be very low\n');
fprintf('Hence, the price of the DIBC is very low\n');