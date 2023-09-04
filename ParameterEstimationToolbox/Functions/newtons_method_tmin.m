function x = newtons_method_tmin(FunJac, xinit, tol, maxit, varargin)

k = 0;
x = xinit;
[f,J] = feval(FunJac,x,varargin{:});
while( (k < maxit) && (f > tol) )
    k = k+1;
    x = x - f/J;
    [f,J] = feval(FunJac,x,varargin{:});
end

