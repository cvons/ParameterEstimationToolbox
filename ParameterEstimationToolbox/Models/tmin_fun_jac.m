function [f,J] = tmin_fun_jac(x,p)

T = x;
V = p.V;
k0 = p.k0;
EaR = p.EaR;
CAin = p.CAin;
CBin = p.CBin;
Tin = p.Tin;
Fmax = p.Fmax;
beta = p.beta;

f = beta*k0*exp(-EaR/T)*(CAin+(Tin-T)/beta)*(CBin+2*(Tin-T)/beta)-((T-Tin)*Fmax)/V;

J = (beta*k0*EaR*exp(-EaR/T)*(CAin+(Tin-T)/beta)*(CBin+2*(Tin-T)/beta))/T^2 ...
    -k0*exp(-EaR/T)*(CBin+2*(Tin-T)/beta)-2*k0*exp(-EaR/T)*(CAin+(Tin-T)/beta)-Fmax/V;