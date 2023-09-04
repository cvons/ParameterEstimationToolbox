function x = newtons_method_ode(FunJac, tk, xk, dt, xinit, tol, maxit, varargin)
% NewtonsMethodODE: Solve for the next state using Newton's method.
%
% Inputs:
%   - FunJac: A function handle representing the ODE and its Jacobian.
%   - tk: Current time.
%   - xk: Current state.
%   - dt: Time step.
%   - xinit: Initial guess for the next state.
%   - tol: Tolerance for convergence.
%   - maxit: Maximum number of iterations.
%   - varargin: Additional arguments passed to FunJac.
%
% Output:
%   - x: The next state after applying Newton's method.

k = 0;            % Initialize iteration count
t = tk + dt;      % Calculate the next time
x = xinit;        % Initialize state variable

% Compute the derivative and Jacobian using FunJac
[f, J] = feval(FunJac, t, x, varargin{:});

% Calculate the residual R
R = x - f * dt - xk;

I = eye(length(xk));  % Identity matrix of appropriate size

% Newton's method iteration
while (k < maxit) && (norm(R, 'inf') > tol)
    k = k + 1;
    
    % Compute the derivative of the residual with respect to x (dR/dx)
    dRdx = I - J * dt;
    
    % Solve for the update in the state variables
    dx = dRdx \ R;
    
    % Update the state
    x = x - dx;
    
    % Recalculate the derivative and Jacobian with the updated state
    [f, J] = feval(FunJac, t, x, varargin{:});
    
    % Recalculate the residual R
    R = x - dt * f - xk;
end