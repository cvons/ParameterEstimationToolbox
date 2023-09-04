function [p_est, Y, X, T, hessian] = param_est_sim(fun_jac_drift, fun_diff, ...
    meas_fun_jac, U, time_interval, N, x0, p, p_guess, obj_fun_type, ...
    x0_est, P0_est, lb, ub, seed)
% param_est_sim Simulates a system, performs measurements, and estimates parameters.
%
% This function simulates a system using the provided drift and diffusion
% functions, performs measurements using the measurement function and
% estimates parameters using measurement data. It returns the estimated
% parameters, measurement data, simulated states and time points.
%
% Inputs:
%   fun_jac_drift: Function handle for the drift function of the system.
%       Should return both the drift function and its Jacobian.
%   fun_diff: Function handle for the diffusion function of the system
%   meas_fun_jac: Function handle for the measurement function and its Jacobian.
%   U: Input trajectory for the system. If the model has no input, provide an empty array, [].
%   time_interval: Time interval [t_start, t_finish] for simulation
%   N: Total number of time steps for simulation
%   x0: Initial state vector for simulation
%   p: Struct of system parameters (p.param1, p.param2, etc.)
%   p_guess: Struct of initial parameter guesses (p.param2, p.param4, etc.)
%   obj_fun_type: Type of objective function for parameter estimation
%   x0_est: Initial state vector for estimation
%   P0_est: Initial covariance matrix for estimation
%   seed: Seed for the random number generator (optional)
%
% Note: Parameters 'p' and 'p_guess' should be provided as a struct, i.e., p.param1, p.param2, etc.
%
% Outputs:
%   p_est: Estimated parameter struct
%   Y: Simulated measurement data
%   X: Simulated state trajectory
%   T: Time points for simulation and measurement
%   hessian: Hessian matrix of the estimation (optional, only returned if requested).
%
        
    % Check if seed argument is provided, otherwise set it to an empty array
    if nargin < 15
        seed = [];
    end

    % Assertions for input types and dimensions
    assert(isa(fun_jac_drift, 'function_handle'), 'fun_jac_drift should be a function handle');
    assert(isa(fun_diff, 'function_handle'), 'fun_diff should be a function handle');
    assert(isa(meas_fun_jac, 'function_handle'), 'meas_fun_jac should be a function handle');
    assert(isnumeric(U) && numel(U) == N, 'U should be a numerical array of length N');
    assert(isnumeric(time_interval) && numel(time_interval) == 2, 'time_interval should be a numerical array of 2 elements');
    assert(isscalar(N) && N > 0 && round(N) == N, 'N should be a positive integer');
    assert(isnumeric(x0) && iscolumn(x0), 'x0 should be a column vector');
    assert(isstruct(p), 'p should be a struct');
    assert(isfield(p, 'R'), 'The struct p must have a field named R.');
    assert(isstruct(p_guess), 'p_guess should be a struct');
    assert(ismember(obj_fun_type, {'ls', 'ml', 'map'}), 'obj_fun_type should be either ''ls'', ''ml'', or ''map''');
    assert(isnumeric(x0_est) && iscolumn(x0_est), 'x0_est should be a column vector');
    assert(isnumeric(P0_est) && ismatrix(P0_est) && size(P0_est, 1) == numel(x0_est) && size(P0_est, 2) == numel(x0_est), ...
        'P0_est should be a square matrix of size numel(x0_est) x numel(x0_est)');
    assert(isnumeric(lb) && isnumeric(ub) && numel(lb) == numel(ub) && numel(lb) == numel(fieldnames(p_guess)), ...
        'lb and ub should be vectors of the same length as the number of parameters that needs to be estimated');

    [Y, X, T] = simulate_model(fun_jac_drift, fun_diff, meas_fun_jac, ...
        U, time_interval, N, x0, p, seed);

    [p_est, hessian] = estimate_parameters(fun_jac_drift, fun_diff, Y, U, x0_est, P0_est, ...
        T, p, p_guess, meas_fun_jac, obj_fun_type, lb, ub);

end