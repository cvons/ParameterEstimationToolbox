function [p_est, hessian] = estimate_parameters(fun_jac_drift, fun_diff, Y, U, x0_est, P0_est, T, ...
    p, p_guess, meas_fun_jac, obj_fun_type, lb, ub)
% estimate_parameters: Estimate model parameters using measured data.
%
% This function estimates model parameters using measured data and an optimization process.
%
% Inputs:
%   - fun_jac_drift: Function handle for the drift function of the system. Should return both the drift function and its Jacobian.
%   - fun_diff: Function handle for the diffusion function of the system.
%   - Y: Measured data.
%   - U: Input trajectory for the system. If the model has no input, provide an empty array, [].
%   - x0_est: Initial guess for the state vector.
%   - P0_est: Initial covariance matrix for the state vector.
%   - T: Time points for measurements.
%   - p: Struct of system parameters.
%   - p_guess: Struct of initial parameter guesses.
%   - meas_fun_jac: Function handle for the measurement function and its Jacobian.
%   - obj_fun_type: Type of objective function for parameter estimation ('ls', 'ml', or 'map').
%   - lb: Lower bounds for estimated parameters.
%   - ub: Upper bounds for estimated parameters.
%
% Outputs:
%   - p_est: Estimated parameter struct.
%   - hessian: Hessian matrix of the estimation (optional, only returned if requested).

    % Assertions for input types and dimensions
    assert(isa(fun_jac_drift, 'function_handle'), 'fun_jac_drift should be a function handle');
    assert(isa(fun_diff, 'function_handle'), 'fun_diff should be a function handle');
    assert(isa(meas_fun_jac, 'function_handle'), 'meas_fun_jac should be a function handle');
    assert(isnumeric(x0_est) && iscolumn(x0_est), 'x0_est should be a column vector');
    assert(isnumeric(P0_est) && ismatrix(P0_est) && size(P0_est, 1) == numel(x0_est) && size(P0_est, 2) == numel(x0_est), ...
        'P0_est should be a square matrix of size numel(x0_est) x numel(x0_est)');
    assert(isstruct(p), 'p should be a struct');
    assert(isfield(p, 'R'), 'The struct p must have a field named R.');
    assert(isstruct(p_guess), 'p_guess should be a struct');
    assert(ismember(obj_fun_type, {'ls', 'ml', 'map'}), 'obj_fun_type should be either ''ls'', ''ml'', or ''map''');
    assert(isnumeric(lb) && isnumeric(ub) && numel(lb) == numel(ub) && numel(lb) == numel(fieldnames(p_guess)), ...
        'lb and ub should be vectors of the same length as the number of parameters that need to be estimated');

    % Create a struct containing only the fields present in p_guess but not in p
    p_fixed = struct_subtract(p, p_guess);

    if nargout == 2
        % Estimate the parameters and their Hessian matrix using the
        % provided functions and inputs
        [p_est, hessian] = parameter_estimation(fun_jac_drift, fun_diff, ...
            Y, U, x0_est, P0_est, T, p_guess, meas_fun_jac, p_fixed, ...
            obj_fun_type, lb, ub);
    else
        % Estimate the parameters using the provided functions and inputs
        p_est = parameter_estimation(fun_jac_drift, fun_diff, Y, U, x0_est, ...
            P0_est, T, p_guess, meas_fun_jac, p_fixed, obj_fun_type, lb, ub);
    end
end