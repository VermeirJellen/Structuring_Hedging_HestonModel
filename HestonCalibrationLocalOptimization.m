function calibratedHestonParameters = HestonCalibrationLocalOptimization(initialGuessHeston,K,T,S0,r,q,optionType,carrMadanPrecomputed,optionMidMarket,weightsOptimization,earlyStopping,maxEval)

kappa=initialGuessHeston(1); eta=initialGuessHeston(2); theta=initialGuessHeston(3);
fellerConstraint=2*kappa*eta-theta^2; % Avoid negative variance: Value must be greater than or equal to zero.
initialGuessFeller = [fellerConstraint initialGuessHeston(2:5)];

lowerBound = [0 0 0 -1 0]; upperBound = [20 1 5 0 1];
if(earlyStopping) % use StopFunction
	optimizationSettings = optimset('OutputFcn',@StopFunction,'TolX',1e-6,'TolFun',1e-6,'MaxFunEvals',10000);
else % do not use StopFunction
	if nargin < 12
		maxEval=10000; % Full Convergence
	end
	optimizationSettings = optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals',maxEval);
end;

[calibratedParams,fval,exitflag,output] = fminsearchbnd(@(paramsToCalibrate) MinimizationCriterionHeston(paramsToCalibrate,K,T,S0,r,q,optionType,carrMadanPrecomputed,optionMidMarket,weightsOptimization),initialGuessFeller,lowerBound,upperBound,optimizationSettings);
calibratedParams(1) = (calibratedParams(1)+calibratedParams(3)^2)/(2*calibratedParams(2)); % extract kappa from fellerConstraint

calibratedHestonParameters = calibratedParams;