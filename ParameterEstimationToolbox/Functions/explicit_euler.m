function [T, X] = explicit_euler(fun, tspan, N, x0, varargin)
% ExplicitEulerFixed implements the Explicit Euler fixed step method.
%
% This function implements the Explicit Euler fixed step method to numerically
% solve a system of ordinary differential equations (ODEs) using the provided
% function handle (fun), time span (tspan), number of steps (N), initial
% state vector (x0) and optional additional arguments (varargin).
% It returns the time vector (T) and state trajectory matrix (X).
%
% Inputs:
%   fun: Function handle for the ODE system
%   tspan: Time span [t0, tN]
%   N: Number of steps
%   x0: Initial state vector
%   varargin: Optional additional arguments for the function
%
% Outputs:
%   T: Time vector
%   X: State trajectory matrix

% Compute step size and allocate memory
t0 = tspan(1);
tN = tspan(2);
dt = (tN - t0) / N;
nx = size(x0, 1);
X = zeros(nx, N + 1);
T = zeros(1, N + 1);

% Explicit Euler Method
T(:, 1) = t0;
X(:, 1) = x0;
for k = 1:N
    f = feval(fun, T(k), X(:, k), varargin{:});
    T(:, k + 1) = T(:, k) + dt;
    X(:, k + 1) = X(:, k) + f * dt;
end