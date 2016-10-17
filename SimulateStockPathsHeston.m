function [stock] = SimulateStockPathsHeston(hestonParameters,S0,T,r,q,nrPeriods,nrPaths,Z1,Z3,maturity)

kappa=hestonParameters(1); 
eta=hestonParameters(2); 
theta=hestonParameters(3); 
corr=hestonParameters(4); 
sig=hestonParameters(5);

M=nrPaths;
N=nrPeriods;
dt = T/N;

if isempty(Z1)
	Z1 = normrnd(0,1,M,N);
end
if isempty(Z3)
	Z3 = normrnd(0,1,M,N);
end

Z2 = corr*Z1 + (sqrt(1-corr^2)*Z3);

zerosMatrix = zeros(M,N); onesVector = ones(M,1);
S = [S0*onesVector, zerosMatrix];
V = [sig^2*onesVector,zerosMatrix];

% Simulation of N-step trajectories for the Variance and price of the underlying asset (Euler Milstein)
for i = 1:N
    V(:,i+1) = V(:,i) + (kappa*(eta-V(:,i))-theta^2/4)*dt + theta*sqrt(V(:,i)).*Z2(:,i)*sqrt(dt) + (theta^2*dt*(Z2(:,i).^2))/4;
    V(:,i+1) = V(:,i+1).*(V(:,i+1)>0); % keep values positive
	S(:,i+1) = S(:,i).*(1+(r-q)*dt + sqrt(V(:,i))*sqrt(dt).*Z1(:,i));
end

if(maturity) % return simulated stock prices at maturity
	stock=S(:,N+1);
else % return full simulated paths
	stock=S;
end