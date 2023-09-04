function [Y, X, T] = simulate_model(fun_jac_drift, fun_diff, ...
    meas_fun, U, time_interval, N, x0, p, seed)
% simulate_model: Simulate a dynamic system and generate measurements.
%
% This function simulates a dynamic system using numerical integration, taking into account
% system drift and diffusion, and then generates measurements based on a measurement function.
%
% Inputs:
%   - fun_jac_drift: Function handle for the drift function of the system. Should return both the drift function and its Jacobian.
%   - fun_diff: Function handle for the diffusion function of the system.
%   - meas_fun: Function handle for the measurement function.
%   - U: Input trajectory for the system. If the model has no input, provide an empty array, [].
%   - time_interval: Time interval [t_start, t_finish] for simulation.
%   - N: Total number of time steps for simulation.
%   - x0: Initial state vector for simulation.
%   - p: Struct of system parameters.
%   - seed: Seed for the random number generator (optional).
%
% Outputs:
%   - Y: Simulated measurement data.
%   - X: Simulated state trajectory.
%   - T: Time points for simulation and measurement.

    % Assertions for input types and dimensions
    assert(isa(fun_jac_drift, 'function_handle'), 'fun_jac_drift should be a function handle');
    assert(isa(fun_diff, 'function_handle'), 'fun_diff should be a function handle');
    assert(isa(meas_fun, 'function_handle'), 'meas_fun should be a function handle');
    assert(isnumeric(U) && numel(U) == N, 'U should be a numerical array of length N');
    assert(isnumeric(time_interval) && numel(time_interval) == 2, 'time_interval should be a numerical array of 2 elements');
    assert(isscalar(N) && N > 0 && round(N) == N, 'N should be a positive integer');
    assert(isnumeric(x0) && iscolumn(x0), 'x0 should be a column vector');
    assert(isstruct(p), 'p should be a struct');
    assert(isfield(p, 'R'), 'The struct p must have a field named R.');

    % Check if seed argument is provided, otherwise set it to an empty array
    if nargin < 9
        seed = [];
    end

    nx = length(x0); % Number of states

    % Simulate the system using the explicit-explicit Euler-Maruyama method
    [W, T, ~] = std_wiener_process(time_interval, N, nx, seed);
    X = euler_maruyama(fun_jac_drift, fun_diff, T, x0, W, U, p);

    % Simulate the measurement using the measurement function
    Y = measurement(meas_fun, X, p.R, seed);
end