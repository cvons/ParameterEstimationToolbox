function X = euler_maruyama(ffun, gfun, T, x0, W, U, varargin)
% EULER_MARUYAMA Simulate a stochastic differential equation using Euler-Maruyama method.
%
% Syntax: X = euler_maruyama(ffun, gfun, T, x0, W, U, varargin)
% X     : State trajectory
%
% ffun  : Function handle for drift function of the SDE
% gfun  : Function handle for diffusion function of the SDE
% T     : Time points
% x0    : Initial state vector
% W     : Wiener process trajectory
% U     : Input trajectory
% varargin : Additional arguments for ffun and gfun

N = length(T); % Number of time steps
nx = length(x0); % Number of states
X = zeros(nx, N); % Initialize state trajectory

X(:, 1) = x0; % Set initial state

if nx == 1 % Scalar case
    for k = 1:N-1
        f = feval(ffun, T(k), X(k), U(k), varargin{:});
        g = feval(gfun, T(k), X(k), U(k), varargin{:});
        dt = T(k+1) - T(k);
        dW = W(k+1) - W(k);
        psi = X(k) + g * dW;
        X(k+1) = psi + f * dt;
    end
else % Multidimensional case
    for k = 1:N-1
        f = feval(ffun, T(k), X(:, k), U(k), varargin{:});
        g = feval(gfun, T(k), X(:, k), U(k), varargin{:});
        dt = T(k+1) - T(k);
        dW = W(:, k+1) - W(:, k);
        psi = X(:, k) + g .* dW;
        X(:, k+1) = psi + f * dt;
    end
end
end
