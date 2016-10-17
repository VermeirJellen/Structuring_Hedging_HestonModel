function stop = StopFunction(x, optimValues, state)
stop = false;
% Check if objective function is less than stopping criterion
global globalStoppingCriterionBidAsk;
if optimValues.fval < globalStoppingCriterionBidAsk
    stop = true;
end