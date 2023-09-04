function xdot = cstr_1d_drift(t,x,u,p)

% Initialize temperature, T, and flow rate, U.
T = x;
F = u;

% Initialize parameters
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

% Define the differential equation
Tdot = (F/V)*(Tin - T) + RT; 
xdot = Tdot;