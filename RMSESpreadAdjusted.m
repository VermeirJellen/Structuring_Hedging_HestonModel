function adjustedRMSE = RMSESpreadAdjusted(modelPrices,marketPrices,bid,ask)

spreadAdjustedError=max(0,max(modelPrices-ask,bid-modelPrices));
adjustedRMSE = sqrt(sum(spreadAdjustedError.^2)/length(marketPrices));