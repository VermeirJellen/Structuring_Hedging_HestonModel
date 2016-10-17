function rmse = RMSE(modelPrices,marketPrices)

rmse = sqrt(sum((modelPrices-marketPrices).^2)/length(marketPrices));