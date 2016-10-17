% Calculate Black Scholes prices and Greeks for a vanilla option or an exotic barrier option, depending on the type input
% Overview of input 'type':
% 0 = call, 2 = down and out barrier call, 3 = down and in barrier call, 4 = up and out barrier call, 5 = up and in barrier call
% 1 = put, 6 = down and out barrier put, 7 = down and in barrier put, 8 = up and out barrier put, 9 = up and in barrier put
function [priceAndGreeks] = BSBarrierPricerAndGreeksAnalytical(sigma,St,K,B,T,r,q,type)

phi = -1; eta = -1;
if(any(type==[0,2,3,4,5]))
	phi=1; % 1 if call
end;
if(any(type==[2,3,6,7])) 
	eta=1; % 1 if barrier approached from above (down and out / down and in)
end;

% Precomputed Variables
tau = T; thetaPlus = (r-q)/sigma + sigma/2; thetaMinus = (r-q)/sigma - sigma/2;
mu = sigma*thetaMinus; lambda= 1+mu/(sigma^2); a=mu/(sigma^2); p=(mu+sigma^2)*tau;
z = 1-K/B; d = exp(-r*tau); f = exp(-q*tau);
X = (log(St/K) + p) / (sigma*sqrt(tau)); % = Black Scholes 'd1'
x1 = (log(St/B)+ p) / (sigma*sqrt(tau));
y = (log(B^2/(St*K)) + p) / (sigma*sqrt(tau));
y1 = ((log(B/St) + p))/(sigma*sqrt(tau));

priceA1 = phi*St*f*normcdf(phi*X) - phi*K*d*normcdf(phi*(X-sigma*sqrt(tau))); % X-sigma*sqrt(tau) = Black-Scholes 'd2'
priceA2 = phi*St*f*normcdf(phi*x1) - phi*K*d*normcdf(phi*(x1-sigma*sqrt(tau)));
priceA3 = (phi*(B/St)^(2*lambda-2))*(St*f*(B/St)^2*normcdf(eta*y) - K*d*normcdf(eta*(y-sigma*sqrt(tau))));
priceA4 = (phi*((-2*mu)/((sigma^2)*St))*((B/St)^(2*lambda-2))) * (St*f*(B/St)^2*normcdf(eta*y1) - K*d*normcdf(eta*(y1-sigma*sqrt(tau)))) - phi*((B/St)^(2*lambda))*f*normcdf(eta*y1) - phi*eta*f*((B/St)^(2*lambda))*normpdf(y1)*z/(sigma*sqrt(tau));

deltaA1 = phi*f*normcdf(phi*X);
deltaA2 = phi*f*normcdf(phi*x1) + f*normpdf(x1)*z/(sigma*sqrt(tau));
deltaA3 = phi*(-2*mu)/(sigma^2*St)*((B/St)^(2*lambda-2)) * (St*f*(B/St)^2*normcdf(eta*y) - K*d*normcdf(eta*(y-sigma*sqrt(tau)))) - phi*((B/St)^(2*lambda))*f*normcdf(eta*y);
deltaA4 = phi*((B/St)^(2*lambda-2)) * (St*f*(B/St)^2*normcdf(eta*y1) - K*d*normcdf(eta*(y1-sigma*sqrt(tau))));

gammaA1 = f*normpdf(X)/(St*sigma*sqrt(tau));
gammaA2 = f*normpdf(x1)/((St*sigma*sqrt(tau))*(1-z*x1/(sigma*sqrt(tau))));
gammaC3 = phi*((B/St)^(2*lambda-2)) * (St*f*((B/St)^2)*normcdf(eta*y) - K*d*normcdf(eta*(y-sigma*sqrt(tau))));
gammaB3 = -2*mu/(sigma^2*St)*gammaC3 - phi*((B/St)^(2*lambda))*f*normcdf(eta*y);
gammaA3 = 2*mu/(sigma^2*St)*(gammaC3/St - gammaB3) + phi*f*(B^(2*lambda))/(St^(2*lambda+1))*(2*lambda*normcdf(eta*y) + eta*normpdf(y)/(sigma*sqrt(tau)));
gammaC4 = phi*((B/St)^(2*lambda-2)) * (St*f*((B/St)^2)*normcdf(eta*y1) - K*d*normcdf(eta*(y1-sigma*sqrt(tau))));
gammaB4 = -2*mu/(sigma^2*St)*gammaC4 - phi*((B/St)^(2*lambda))*f*normcdf(eta*y1) - phi*eta*f*((B/St)^(2*lambda))*normpdf(y1)*z/(sigma*sqrt(tau));
gammaA4 = 2*mu/(sigma^2*St)*(gammaC4/St - gammaB4) + phi*f*(B^(2*lambda))/(St^(2*lambda+1))*(2*lambda*normcdf(eta*y1) + eta*normpdf(y1)/(sigma*sqrt(tau))) + phi*eta*f*z*normpdf(y1)*((B/St)^(2*lambda))/(St*sigma*sqrt(tau))*(2*lambda-y1/(sigma*sqrt(tau)));

thetaA1 = -1/2*sigma*St*f*normpdf(X)/sqrt(tau) + phi*St*f*normcdf(phi*X)*q - phi*K*d*normcdf(phi*(X-sigma*sqrt(tau)))*r;
thetaA2 = -1/2*sigma*St*f*normpdf(x1)*K/(B*sqrt(tau)) + phi*St*f*normcdf(phi*x1)*q - phi*K*d*normcdf(phi*(x1-sigma*sqrt(tau)))*r - St*f*normpdf(x1)*z*y1/(2*tau);
thetaA3 = -phi*((B/St)^(2*lambda))*St*f*eta*normpdf(y)*(1/2)*sigma/sqrt(tau) + phi*((B/St)^(2*lambda-2))*(q*St*f*((B/St)^2)*normcdf(eta*y) - r*K*d*normcdf(eta*(y-sigma*sqrt(tau))));
thetaA4 = -phi*((B/St)^(2*lambda))*St*f*eta*normpdf(y1)*(x1/(2*tau)*z + (1/2)*sigma*K/(sqrt(tau)*B)) + phi*((B/St)^(2*lambda-2))*(q*St*f*((B/St)^2)*normcdf(eta*y1) - r*K*d*normcdf(eta*(y1-sigma*sqrt(tau))));

vegaA1 = St*f*normpdf(X)*sqrt(tau);
vegaA2 = St*f*normpdf(x1)*(sqrt(tau)-x1*z/sigma);
vegaB3 = phi*((B/St)^(2*lambda-2))*(St*f*((B/St)^2)*normcdf(eta*y) - K*d*normcdf(eta*(y-sigma*sqrt(tau))));
vegaA3 = -4/(sigma^3)*log(B/St)*(r-q)*vegaB3 + phi*((B/St)^(2*lambda))*St*f*eta*normpdf(y)*sqrt(tau);
vegaB4 = phi*((B/St)^(2*lambda-2))*(St*f*((B/St)^2)*normcdf(eta*y1) - K*d*normcdf(eta*(y1-sigma*sqrt(tau))));
vegaA4 = -4/(sigma^3)*log(B/St)*(r-q)*vegaB4 + phi*((B/St)^(2*lambda))*St*f*eta*normpdf(y1)*((sqrt(tau)-y1/sigma)*z + K/B*sqrt(tau));

if (type==0 | type==1 | (type==5 & K > B) | (type==7 & K <= B))
	priceAndGreeks = [priceA1, deltaA1, gammaA1, thetaA1, vegaA1];
elseif ((type==5 & K<=B) | (type == 7 & K>B))
	priceAndGreeks = [priceA2-priceA3+priceA4, deltaA2-deltaA3+deltaA4, gammaA2-gammaA3+gammaA4, thetaA2-thetaA3+thetaA4, vegaA2-vegaA3+vegaA4];
elseif ((type==9 & K>B) | (type==3 & K<=B))
	priceAndGreeks = [priceA1-priceA2+priceA4, deltaA1-deltaA2+deltaA4, gammaA1-gammaA2+gammaA4, thetaA1-thetaA2+thetaA4, vegaA1-vegaA2+vegaA4];
elseif ((type==9 & K<=B) | (type==3 & K > B))
	priceAndGreeks = [priceA3, deltaA3, gammaA3, thetaA3, vegaA3];
elseif ((type== 4 & K>B) | (type==6 & K<=B))
	priceAndGreeks = [0,0,0,0,0];
elseif ((type==4 & K<=B) | (type==6 & K>B))
	priceAndGreeks = [priceA1-priceA2+priceA3-priceA4, deltaA1-deltaA2+deltaA3-deltaA4, gammaA1-gammaA2+gammaA3-gammaA4, thetaA1-thetaA2+thetaA3-thetaA4, vegaA1-vegaA2+vegaA3-vegaA4];
elseif ((type == 8 & K>B) | (type==2 & K<=B))
	priceAndGreeks = [priceA2-priceA4, deltaA2-deltaA4, gammaA2-gammaA4, thetaA2-thetaA4, vegaA2-vegaA4];
elseif ((type == 8 & K<=B) | (type==2 & K>B))
	priceAndGreeks = [priceA1-priceA3, deltaA1-deltaA3, gammaA1-gammaA3, thetaA1-thetaA3, vegaA1-vegaA3];
end;

% 0 = call
% 1 = put

% 2 = down and out barrier call
% 3 = down and in barrier call
% 4 = up and out barrier call
% 5 = up and in barrier call

% 6 = down and out barrier put
% 7 = down and in barrier put
% 8 = up and out barrier put
% 9 = up and in barrier put