clc; clear;

fprintf('Reading market data..\n');

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

fprintf('Starting Heston calibration process..\n');
%% HESTON MODEL CALIBRATION (CARR MADAN WITH CROSS VALIDATION APPROACH %%
% Carr Madan Parameters
N=4096; alpha=1.5; gridSpace=0.25; simpsonIntegrand=1;
carrMadanPrecomputed = PrecomputationCarrMadanParameters(N,alpha,gridSpace,simpsonIntegrand); % Precompute modelindependent carr-madan parameters and variables

% Initial educated guess for HestonParameters
kappa=7; eta=0.7; theta=0.7; corr=-0.5; sig=0.18;
initialGuessHeston = [kappa eta theta corr sig];

% 3. Use Implied volatility as weights. Use out of sample cross validation for calibration.
weightsOptimization = 1./BlackScholesImpliedVolatility(S0,K,r,q,T,optionType,optionMidMarket,false);
% weightsOptimization = 1./abs(ask-bid);
% weightsOptimization = ones(1,length(K));
[hestonParameters, outOfSampleRMSE, outOfSampleRMSEAdjusted] = CalibrateHestonParametersCrossValidation(marketData,carrMadanPrecomputed,initialGuessHeston,weightsOptimization,true); % (14.3 seconds)
% outOfSampleRMSE1 = 0.1367 % outOfSampleAdjustedRMSE1 = 0.0163

hestonParameters;
outOfSampleRMSE;
outOfSampleRMSEAdjusted;

fprintf('Calibrated HestonParameters:\n');
fprintf('kappa: %.4f \n',hestonParameters(1));
fprintf('eta: %.4f \n',hestonParameters(2));
fprintf('theta: %.4f \n',hestonParameters(3));
fprintf('rho: %.4f \n',hestonParameters(4));
fprintf('sig0: %.4f \n\n',hestonParameters(5));

barrierLevel=0.70; KMaturity=S0; nrDays = 652; % 652 days until maturity
rMaturity=spline(interestRates(:,1), interestRates(:,2), nrDays); TMaturity=nrDays/365;
confidence=0.05; % parameter for confidence intervals

% Control variate AND anti thetic variable, 75000 pairs
nrPaths=100000; nrPeriods = nrDays*2; % 
tic; [priceDOBCAntiTheticControlVariate, lowerLevel, upperLevel] = MonteCarloPricerDOBCHeston(hestonParameters,KMaturity,TMaturity,S0,barrierLevel,rMaturity,q,nrPeriods,nrPaths,confidence,true,true); toc; % 3.6 seconds
confidenceInterval = upperLevel-lowerLevel; % confidence interval @ 95 percent (upperLevel-lowerLevel) =  0.0206
priceDOBCAntiTheticControlVariate; %  We use 4.50, small variations each run

fprintf('Monte Carlo:\n');
fprintf('Price: %.4f\n',priceDOBCAntiTheticControlVariate);
fprintf('95 percent confidence interval: %.4f\n',confidenceInterval);

hestonPriceDOBC = 4.50;
fprintf('In our calculations we use %.4f as our DOBC price\n\n',hestonPriceDOBC);
customerPriceDOBC = hestonPriceDOBC; % price at upper level in calculations. Take uncertainty of the exotic price/hedging into account.


%%% PPPN CALCULATIONS %%%
N = 1000000;
principalProtection = 0.95;
bankAccount = N*principalProtection*exp(-rMaturity*(nrDays/365));
remainingNotionalAfterProtection = N-bankAccount;

fixedProfit = N*0.01; % Profit Margin
remainingNotionalAfterProfitTaking = remainingNotionalAfterProtection-fixedProfit;

nrDOBCSold = floor(remainingNotionalAfterProfitTaking/customerPriceDOBC);

% Participation "on the upside" is the amount of options we buy versus the amount of stocks we could have bought with the original notion;
notionalStocks = N/S0;
participation = nrDOBCSold/notionalStocks;

fprintf('Structuring process:\n');
fprintf('Notional is 1000000USD\n');
fprintf('We put %.2fUSD on the bankaccount',bankAccount);
fprintf('We take %.2fUSD as our profitmargin\n',fixedProfit);
fprintf('We now have %.2fUSD left to buy DOBC options\n',remainingNotionalAfterProfitTaking);
fprintf('From this amount we can buy %d DOBC options\n',nrDOBCSold);
fprintf('This corrsponds to a participation rate of %.4f\n',participation);
fprintf('Investor pnl at maturity is plotted on the graph\n\n');

% Plot payoff 
pnlStock = @(St) N/S0*St-N;
pnlRiskFree= N*exp(rMaturity*(nrDays/365))-N;

payoffBankAccount = bankAccount*exp(rMaturity*(nrDays/365));
% Note: barrier not breached when ST > 0. Obviously breached when ST < barrier.
%pnlPPPBarrierNotBreached = @(St) payoffBankAccount + participation*N*max((St-S0)/S0,0) - N;
pnlPPPBarrierNotBreached = @(St) payoffBankAccount + max(nrDOBCSold*(St-S0),0) - N; % Strike of DOBC = S0;
pnlPPPBarrierBreached = @(St) payoffBankAccount - N;

%%%%%%%%%%%%%%%%%%%%%%% Graph the results and perform some analysis *****************************
St=0:0.001:100; % Stock price
intersectIndex = find(pnlPPPBarrierNotBreached(St)-pnlStock(St) < eps,1); % Find pnl intersection (pnlPPP - pnlStock == 0)
intersectX = St(intersectIndex); intersectY = pnlPPPBarrierNotBreached(intersectX); % X and Y coordinates of intersection
breakEvenIndex = find(abs(pnlPPPBarrierNotBreached(St)) < 50,1); % find breakevenindex
breakEvenX = St(breakEvenIndex); % This is the stock price St at maturity for which the holder of the PPP will break even
% intersectRiskFreeIndex = find(pnlPPPBarrierNotBreached(St)-pnlRiskFree > eps,1);
% intersectRiskFree = St(intersectRiskFreeIndex); % This is the stock price St at maturity for which the holder of the RC will generate the same profit as the risk free rate (bank)

figure;
hold on;
h0 = graph2d.constantline(payoffBankAccount-N,'LineStyle','--','Color','b');
changedependvar(h0,'y');

hPnl = plot(St,pnlPPPBarrierNotBreached(St),'-',St,pnlStock(St),'-');
title(sprintf('Investment PnL at maturity (N = %.0fUSD)',N));
xlabel('Stock price St at time T (USD)') % x-axis label
ylabel('Pnl at maturity (USD)') % y-axis label

%hRiskFree = graph2d.constantline(pnlRiskFree, 'LineStyle','-', 'Color','r'); %changedependvar(hRiskFree,'y');
%hBarrier = graph2d.constantline(B*S0, 'LineStyle','--', 'Color','r'); %changedependvar(hBarrier,'x');
hS0 = graph2d.constantline(S0, 'LineStyle','--', 'Color','g');
changedependvar(hS0,'x');
hPnlIntersect = plot(intersectX,intersectY,'ko','MarkerSize',8);
hBreakeven = plot(breakEvenX,0,'ro','MarkerSize',8);

legend(sprintf('PPP (Barrier breached)'),'PPP (Barrier not breached)','XLB index',sprintf('S0 = %.2fUSD',S0),sprintf('St = %.2fUSD, pnl = %.2fUSD',intersectX,intersectY),sprintf('St = %.2fUSD, pnl = %.2fUSD',breakEvenX,0),'Location','northwest');
h0 = graph2d.constantline(0,'LineStyle','-', 'Color','r');
changedependvar(h0,'y');
%%%%%%%%%%%%%%%%%%%%%%% Graph the results and perform some analysis *****************************



fprintf('Hedging of the DOBC\n');

%%%%%%%%%%%%%%%%%%%%%%%% HEDGING AT TIME 0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BlackScholes = BSBarrierPricer(sigma,St,K,B,T,r,q,type)
sigmaCalibrated = BlackScholesImpliedVolatilityExotic(S0,KMaturity,barrierLevel*S0,TMaturity,rMaturity,q,2,hestonPriceDOBC); % Implied volatility for DOBC matching heston price
BSPriceAndGreeksDOBC = BSBarrierPricerAndGreeksAnalytical(sigmaCalibrated,S0,KMaturity,barrierLevel*S0,TMaturity,rMaturity,q,2);

fprintf('The black scholes implied volatility for our exotic option is %.2f\n',sigmaCalibrated);
fprintf('The greeks of the DOBC under the BS model are as follows\n');
fprintf('Delta = %.4f\n',BSPriceAndGreeksDOBC(2));
fprintf('Gamma = %.4f\n',BSPriceAndGreeksDOBC(3));
fprintf('Theta = %.4f\n',BSPriceAndGreeksDOBC(4));
fprintf('Vega = %.4f\n\n',BSPriceAndGreeksDOBC(5));

soldOptionPriceAndGreeks = nrDOBCSold*BSPriceAndGreeksDOBC;
fprintf('Total exposure to hedge for shorting the PPPN containing %d DOBC options:\n',nrDOBCSold);
fprintf('Delta = %.4f\n',soldOptionPriceAndGreeks(2));
fprintf('Gamma = %.4f\n',soldOptionPriceAndGreeks(3));
fprintf('Theta = %.4f\n',soldOptionPriceAndGreeks(4));
fprintf('Vega = %.4f\n\n',soldOptionPriceAndGreeks(5));

%% HEDGING INSTRUMENTS %%
%% DELTA VEGA HEDGING and DELTA-VEGA-GAMMA-THETA hedging %%
nrO=146; % (OTM put, K=48.5, T=38 days)
impliedVolHedgingInstrument1 = BlackScholesImpliedVolatility(S0,K(nrO),r(nrO),q,T(nrO),optionType(nrO),optionMidMarket(nrO),false);
hedgingInstrument1 = BSBarrierPricerAndGreeksAnalytical(impliedVolHedgingInstrument1,S0,K(nrO),S0,T(nrO),r(nrO),q,optionType(nrO));
nrO=231; % (ATM call, K=51, T=241 days)
impliedVolHedgingInstrument2 = BlackScholesImpliedVolatility(S0,K(nrO),r(nrO),q,T(nrO),optionType(nrO),optionMidMarket(nrO),false);
hedgingInstrument2 = BSBarrierPricerAndGreeksAnalytical(impliedVolHedgingInstrument2,S0,K(nrO),S0,T(nrO),r(nrO),q,optionType(nrO));
nrO=289; % (OTM call, K=65, T=612 days)
impliedVolHedgingInstrument3 = BlackScholesImpliedVolatility(S0,K(nrO),r(nrO),q,T(nrO),optionType(nrO),optionMidMarket(nrO),false);
hedgingInstrument3 = BSBarrierPricerAndGreeksAnalytical(impliedVolHedgingInstrument3,S0,K(nrO),S0,T(nrO),r(nrO),q,optionType(nrO));
nrO=179; % (OTM put, K=45, T=122 days)
impliedVolHedgingInstrument4 = BlackScholesImpliedVolatility(S0,K(nrO),r(nrO),q,T(nrO),optionType(nrO),optionMidMarket(nrO),false);
hedgingInstrument4 = BSBarrierPricerAndGreeksAnalytical(impliedVolHedgingInstrument4,S0,K(nrO),S0,T(nrO),r(nrO),q,optionType(nrO));
hedgingInstrument5 = [S0, 1, 0, 0, 0]; % The stock

fprintf('We consider following hedging instruments:\n');
fprintf('OTM put, K=48.5, T=38 days. price=%.4f,Delta=%.4f,Gamma=%.4f,Theta=%.4f,Vega=%.4f\n',hedgingInstrument1(1),hedgingInstrument1(2),hedgingInstrument1(3),hedgingInstrument1(4),hedgingInstrument1(5));
fprintf('ATM call, K=51, T=241 days. price=%.4f,Delta=%.4f,Gamma=%.4f,Theta=%.4f,Vega=%.4f\n',hedgingInstrument2(1),hedgingInstrument2(2),hedgingInstrument2(3),hedgingInstrument2(4),hedgingInstrument2(5));
fprintf('OTM call, K=65, T=612 days. price=%.4f,Delta=%.4f,Gamma=%.4f,Theta=%.4f,Vega=%.4f\n',hedgingInstrument3(1),hedgingInstrument3(2),hedgingInstrument3(3),hedgingInstrument3(4),hedgingInstrument3(5));
fprintf('OTM put, K=45, T=122 days. price=%.4f,Delta=%.4f,Gamma=%.4f,Theta=%.4f,Vega=%.4f\n',hedgingInstrument4(1),hedgingInstrument4(2),hedgingInstrument4(3),hedgingInstrument4(4),hedgingInstrument4(5));
fprintf('Stock. price=%.4f,Delta=%.4f,Gamma=%.4f,Theta=%.4f,Vega=%.4f\n\n',hedgingInstrument5(1),hedgingInstrument5(2),hedgingInstrument5(3),hedgingInstrument5(4),hedgingInstrument5(5));

%% DELTA HEDGING WITH STOCK at time 0%%
amountSold = N-fixedProfit;
DOBCPosition = -soldOptionPriceAndGreeks(1); % We sell the option to the client at time 0
nrStocks = floor(soldOptionPriceAndGreeks(2));
stockPosition = nrStocks*S0; % We have to buy delta stock to hedge the DOBC
bankPosition = -(DOBCPosition+stockPosition); % We have to borrow money to buy the stock
totalPositionDeltaHedging = DOBCPosition+stockPosition+bankPosition;

deltaHedgingInstruments = [hedgingInstrument5([2])'];
quantitiesHedgingInstruments = deltaHedgingInstruments\soldOptionPriceAndGreeks(2)';
DOBCTotalValueAndGreeks = -soldOptionPriceAndGreeks;
hedgingPortfolioTotalValueAndGreeks = ([hedgingInstrument5']*quantitiesHedgingInstruments)';
combinationPortfolioTotalValueAndGreeks = DOBCTotalValueAndGreeks + hedgingPortfolioTotalValueAndGreeks;

bankAccount = -combinationPortfolioTotalValueAndGreeks(1); % The total value of the portfolio is negative, which means we have money left to put on the bankaccount
thetaBankAccount = bankAccount*rMaturity; % Very small theta
combinationPortfolioTotalValueAndGreeks(1) = combinationPortfolioTotalValueAndGreeks(1) + bankAccount;
combinationPortfolioTotalValueAndGreeks(4) = combinationPortfolioTotalValueAndGreeks(4) + thetaBankAccount;
finalPortfolio = combinationPortfolioTotalValueAndGreeks;

fprintf('DELTA HEDGIN£G WITH STOCK:\n');
fprintf('We buy %.2f stocks\n',quantitiesHedgingInstruments(1));
fprintf('We borrow an additional %.2fUSD in order to finance the stock purchase\n',-bankAccount);
fprintf('Total exposure for the combined portfolio:\n');
fprintf('pnl = %.4f\n',finalPortfolio(1));
fprintf('Delta = %.4f\n',finalPortfolio(2));
fprintf('Gamma = %.4f\n',finalPortfolio(3));
fprintf('Theta = %.4f\n',finalPortfolio(4));
fprintf('Vega = %.4f\n',finalPortfolio(5));
fprintf('View report for more details\n\n');

%% DELTA VEGA HEDGING (instrument 2 en 5 + bankaccount for value matching) %%
deltaVegaHedgingInstruments = [hedgingInstrument2([2,5])', hedgingInstrument5([2,5])'];
quantitiesHedgingInstruments = deltaVegaHedgingInstruments\soldOptionPriceAndGreeks([2,5])';
DOBCTotalValueAndGreeks = -soldOptionPriceAndGreeks;
hedgingPortfolioTotalValueAndGreeks = ([hedgingInstrument2', hedgingInstrument5']*quantitiesHedgingInstruments)';
combinationPortfolioTotalValueAndGreeks = DOBCTotalValueAndGreeks + hedgingPortfolioTotalValueAndGreeks;

bankAccount = -combinationPortfolioTotalValueAndGreeks(1); % The total value of the portfolio is negative, which means we have money left to put on the bankaccount
thetaBankAccount = bankAccount*rMaturity; % Very small theta
combinationPortfolioTotalValueAndGreeks(1) = combinationPortfolioTotalValueAndGreeks(1) + bankAccount;
combinationPortfolioTotalValueAndGreeks(4) = combinationPortfolioTotalValueAndGreeks(4) + thetaBankAccount;
finalPortfolio = combinationPortfolioTotalValueAndGreeks;

fprintf('DELTA VEGA HEDGIN£G WITH ATM CALL AND STOCK:\n');
fprintf('We buy %.2f stocks\n',quantitiesHedgingInstruments(2));
fprintf('We buy %.2f ATM call options\n',quantitiesHedgingInstruments(1));
fprintf('We have %.2fUSD left to put on the bankaccount\n',bankAccount);
fprintf('Total exposure for the combined portfolio:\n');
fprintf('pnl = %.4f\n',finalPortfolio(1));
fprintf('Delta = %.4f\n',finalPortfolio(2));
fprintf('Gamma = %.4f\n',finalPortfolio(3));
fprintf('Theta = %.4f\n',finalPortfolio(4));
fprintf('Vega = %.4f\n',finalPortfolio(5));
fprintf('View report for more details\n\n');

%% DELTA VEGA IN COMBINATION PORTFOLIO VALUE MATHING %%
deltaVegaHedgingInstruments = [hedgingInstrument1([1,2,5])', hedgingInstrument3([1,2,5])', hedgingInstrument5([1,2,5])'];
quantitiesHedgingInstruments = deltaVegaHedgingInstruments\soldOptionPriceAndGreeks([1,2,5])';
DOBCTotalValueAndGreeks = -soldOptionPriceAndGreeks;
hedgingPortfolioTotalValueAndGreeks = ([hedgingInstrument1', hedgingInstrument3', hedgingInstrument5']*quantitiesHedgingInstruments)';
finalPortfolio = DOBCTotalValueAndGreeks + hedgingPortfolioTotalValueAndGreeks; % 0 portfolio value, delta and vega neutral. gamma short, theta long.

fprintf('DELTA VEGA HEDGIN£G WITH OTM CALL, OTM PUT AND STOCK (Match DOBC portfolio value):\n');
fprintf('We buy %.2f stocks\n',quantitiesHedgingInstruments(3));
fprintf('We buy %.2f OTM put options\n',quantitiesHedgingInstruments(1));
fprintf('We buy %.2f OTM call options\n',quantitiesHedgingInstruments(2));
fprintf('Total exposure for the combined portfolio:\n');
fprintf('pnl = %.4f\n',finalPortfolio(1));
fprintf('Delta = %.4f\n',finalPortfolio(2));
fprintf('Gamma = %.4f\n',finalPortfolio(3));
fprintf('Theta = %.4f\n',finalPortfolio(4));
fprintf('Vega = %.4f\n',finalPortfolio(5));
fprintf('View report for more details\n\n');

%% DELTA-VEGA-GAMMA-THETA hedging %%
priceDeltaGammaThetaVegaHedgingInstruments = [hedgingInstrument1', hedgingInstrument2', hedgingInstrument3', hedgingInstrument4' hedgingInstrument5']; %price 
mDeterminant = det(priceDeltaGammaThetaVegaHedgingInstruments); % non zero
% Note that we also include the actual portfolio value in the linear system: This ensures that hedging portfolio value matches the DOBC portfolio value
% (This Avoids extra theta exposure when borrowing money or adding money to the bankaccount in the scenario where the portfolio values do not match)
quantitiesHedgingInstruments = priceDeltaGammaThetaVegaHedgingInstruments\soldOptionPriceAndGreeks(1:5)';

% CHECK RESULTS %
nrDOBCOriginalPortfolio = -nrDOBCSold; % short the DOBC
nrInstrumentsHedgingPortfolio = quantitiesHedgingInstruments'; % long the hedging portfolio
DOBCTotalValueAndGreeks = (-soldOptionPriceAndGreeks);
HedgingPortfolioTotalValueAndGreeks = (priceDeltaGammaThetaVegaHedgingInstruments*quantitiesHedgingInstruments)';
finalPortfolio = DOBCTotalValueAndGreeks + HedgingPortfolioTotalValueAndGreeks; % zero with small rounding error

fprintf('DELTA VEGA GAMMA THETA HEDGIN£G WITH ALL 4 HEDGING OPTIONS AND STOCK:\n');
fprintf('We buy %.2f stocks\n',quantitiesHedgingInstruments(5));
fprintf('We buy %.2f of the first type of OTM put options\n',quantitiesHedgingInstruments(1));
fprintf('We buy %.2f ATM call options\n',quantitiesHedgingInstruments(2));
fprintf('We buy %.2f OTM call options\n',quantitiesHedgingInstruments(3));
fprintf('We buy %.2f of the second type of OTM put options\n',quantitiesHedgingInstruments(4));
fprintf('Total exposure for the combined portfolio:\n');
fprintf('pnl = %.4f\n',finalPortfolio(1));
fprintf('Delta = %.4f\n',finalPortfolio(2));
fprintf('Gamma = %.4f\n',finalPortfolio(3));
fprintf('Theta = %.4f\n',finalPortfolio(4));
fprintf('Vega = %.4f\n',finalPortfolio(5));
fprintf('View report for more details\n\n');