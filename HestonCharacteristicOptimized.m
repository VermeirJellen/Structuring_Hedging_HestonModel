function [y] = HestonCharacteristicOptimed(u,t,r,q,p)
% p = struct('x0',log(S0),'d',d,'g',g,'param2',param2,'param3',param3,'param4',param4);

p1 = exp(-p.d.*t); % evaluate exp(-dt) only one time
p2 = 1-p.g.*p1;

A = 1i*u*(p.x0+(r-q).*t);
B = p.param3*(p.param2.*t-2*log(p2./(1-p.g)));
C = p.param4*p.param2.*(1-p1)./p2;
y = exp(A+B+C);

end