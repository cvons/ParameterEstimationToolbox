function Tmin = compute_tmin(p)

tol = 1.0e-12;
maxit = 1000;

TminInit = (p.V*p.beta*p.k0*exp(-p.EaR/p.Tin)*p.CAin*p.CBin)/p.Fmax+p.Tin;

Tmin = NewtonsMethod(@fTminFunJac,TminInit,tol,maxit,p);