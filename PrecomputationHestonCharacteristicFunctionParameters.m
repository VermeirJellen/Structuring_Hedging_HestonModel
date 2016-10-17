% Precompute characteristic function related variables that are independent of a particular maturity
% These parameters can be precomputed only 1 time when pricing sets of time dependent options, for performance reasons.
function characteristicFunctionPrecomputedVariables = PrecomputationHestonCharacteristicFunctionParameters(hestonParameters,u,S0);

kappa=hestonParameters(1);
theta=hestonParameters(3);

param1 = kappa-hestonParameters(4)*theta*u*1i; d = sqrt(param1.^2 - theta^2*(-1i*u-u.^2));
param2 = param1-d;
g = param2./(param1+d);
param3 = kappa*hestonParameters(2)/(theta^2);
param4 = hestonParameters(5)^2/theta^2;

characteristicFunctionPrecomputedVariables = struct('x0',log(S0), ...
													'd',d, ...
													'g',g, ...
													'param2',param2, ...
													'param3',param3, ...
													'param4',param4);