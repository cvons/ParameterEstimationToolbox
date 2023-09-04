function [T, X] = implicit_euler(funJac, ta, tb, N, xa, varargin)
% ImplicitEulerFixed: Solve an ODE using the Implicit Euler method with a fixed step size.
%
% Inputs:
%   - funJac: A function handle representing the ODE and its Jacobian.
%   - ta: Initial time.
%   - tb: Final time.
%   - N: Number of time steps.
%   - xa: Initial condition.
%   - varargin: Additional arguments passed to funJac.
%
% Outputs:
%   - T: Array containing the time values at each step.
%   - X: Array containing the state variables at each step.

% Compute step size and allocate memory
dt = (tb - ta) / N;    % Calculate the time step size
nx = size(xa, 1);      % Number of state variables
X = zeros(nx, N + 1);  % Initialize state variable array
T = zeros(1, N + 1);   % Initialize time array

tol = 1.0e-8;           % Tolerance for Newton's method
maxit = 100;            % Maximum number of iterations for Newton's method

% Eulers Implicit Method
T(:, 1) = ta;  % Initial time
X(:, 1) = xa;  % Initial condition
for k = 1:N
    % Compute the derivative using the function funJac
    f = feval(funJac, T(k), X(:, k), varargin{:});
    
    % Update time and state variables using Implicit Euler
    T(:, k + 1) = T(:, k) + dt;
    xinit = X(:, k) + f * dt;
    
    % Solve for the next state using Newton's method
    X(:, k + 1) = newtons_method_ode(funJac, T(:, k), X(:, k), dt, xinit, tol, maxit, varargin{:});
end