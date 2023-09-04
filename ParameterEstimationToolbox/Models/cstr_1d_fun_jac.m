function [f, J] = cstr_1d_fun_jac(t,x,u,p)

% x = T
% u = F
% p = [CAin; CBin; Tin;  ....] i.e all the parameters

T = x;
F = u;

V = p.V;
logk0 = p.logk0;
EaR = p.EaR;
CAin = p.CAin;
CBin = p.CBin;
Tin = p.Tin;

beta = p.beta;

CA = CAin+1/beta*(Tin-T); 
CB = CBin+2/beta*(Tin-T); 

k = exp(logk0)*exp(-EaR/T); 
r = k*CA*CB; 
RT = beta*r; 

Tdot = (F/V)*(Tin - T) + RT; 
f = Tdot;
if nargout == 2
    J = - (F/V) + (beta*k*EaR*CA*CB)/T^2 - k * CB - 2 * k * CA;
end