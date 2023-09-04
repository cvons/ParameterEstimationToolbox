function g = cstr_3d_diffusion(t,x,u,p)

F = u;
V = p.V;
sigma = p.sigma_T;

g = [0;0;F/V*sigma];