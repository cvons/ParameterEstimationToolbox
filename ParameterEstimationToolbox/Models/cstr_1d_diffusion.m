function g = cstr_1d_diffusion(t,x,u,p)

F = u;
V = p.V;
sigma = p.sigma_T;

g = F/V*sigma;