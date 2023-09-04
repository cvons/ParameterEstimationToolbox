function [f, J] = cstr_3d_fun_jac(t,x,u,p)

% x = T
% u = F
% p = [CAin; CBin; Tin;  ....] i.e all the parameters

CA = x(1);
CB = x(2);
T = x(3);
F = u;

V = p.V;
logk0 = p.logk0;
EaR = p.EaR;
CAin = p.CAin;
CBin = p.CBin;
Tin = p.Tin;
beta = p.beta;

k = exp(logk0)*exp(-EaR/T);
r = k*CA*CB; 
RA = -r;
RB = -2*r;
RT = beta*r;

CAdot = (F/V)*(CAin-CA) + RA;
CBdot = (F/V)*(CBin-CB) + RB;
Tdot = (F/V)*(Tin - T) + RT;
f = [CAdot;CBdot;Tdot];

if nargout == 2
    % Calculate the partial derivatives
    dRAdCA = -k * CB;
    dRAdCB = -k * CA;
    dRAdT = -k * EaR/T^2 * CA * CB;
    
    dRBdCA = -2 * k * CB;
    dRBdCB = -2 * k * CA;
    dRBdT = -2 * k * EaR/T^2 * CA * CB;
    
    dRTdCA = beta * -dRAdCA;
    dRTdCB = beta * -dRAdCB;
    dRTdT = (beta*k*EaR*CA*CB)/T^2;

    % Build the Jacobian matrix
    J = [-(F/V) + dRAdCA, dRAdCB, dRAdT;
         dRBdCA, -(F/V) + dRBdCB, dRBdT;
         dRTdCA, dRTdCB, -(F/V) + dRTdT];
end