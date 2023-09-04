function [theta_est, hessian] = parameter_estimation(fun_jac, fun_diff, Y, U, x0, ...
    P0, T, p_guess, meas_fun_jac, p_fixed, obj_type, lb, ub)
% parameter_estimation Estimates parameters using optimization techniques.
%
% This function estimates parameters by optimizing an objective function
% based on the input obj_type using the provided drift function and its Jacobian (fun_jac),
% measurement data (Y), input trajectory (U), initial state (x0),
% prior covariance (P0), time points (T), initial parameter guess struct
% (p_guess), measurement function handle including its Jacobian (meas_fun_jac),
% fixed parameters (p_fixed), optimization objective type (obj_type) and parameter bounds (lb and ub).
% It returns the estimated parameter vector theta_est.
%
% Inputs:
%   fun_jac: Drift function handle including its Jacobian
%   Y: Measurement data
%   U: Input trajectory
%   x0: Initial state
%   P0: Initial covariance
%   T: Time points
%   p_guess: Initial parameter guess struct
%   meas_fun_jac: Measurement function handle including its Jacobian
%   p_fixed: Struct of fixed parameters
%   obj_type: Optimization objective type ('ls', 'ml', or 'map')
%   lb: Lower parameter bounds
%   ub: Upper parameter bounds
%
% Output:
%   theta_est: Estimated parameter vector
    
    % Compute initial state for CD-EKF
    x0 = mvnrnd(x0, P0)';

    % Convert the initial parameter guess struct to an array
    p_guess_fields = fieldnames(p_guess);
    p_initial_vec = zeros(numel(p_guess_fields), 1);
    for i = 1:numel(p_guess_fields)
        p_initial_vec(i) = p_guess.(p_guess_fields{i});
    end

    % Define the objective function based on the input obj_type
    switch obj_type
        case 'ls'
            % Define the least squares objective function
            objective_function = @(theta) ls_objective_function(fun_jac, ...
                fun_diff, theta, Y, U, x0, P0, T, meas_fun_jac, p_fixed, ...
                fieldnames(p_guess));
        case 'ml'
            % Define the maximum likelihood objective function
            objective_function = @(theta) ml_objective_function(fun_jac, ...
                fun_diff, theta, Y, U, x0, P0, T, meas_fun_jac, p_fixed, ...
                fieldnames(p_guess));
        case 'map'
            % Create parameter struct from p_guess and p_fixed
            p = struct_copy(p_guess, p_fixed);

            % Compute P_theta_0
            [~, P_theta_0, ~, ~] = cd_ekf(fun_jac, fun_diff, Y, U, x0, P0, T, meas_fun_jac, p);

            % Provide the prior mean and covariance matrix for MAP estimation
            objective_function = @(theta) map_objective_function(fun_jac, ...
                fun_diff, theta, Y, U, x0, P0, T, meas_fun_jac, p_fixed, ...
                p_guess, P_theta_0(end), fieldnames(p_guess));
    end

    % Set options for the optimization algorithm
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', ...
        'Display', 'iter');

    if nargout == 2
        % Perform optimization to find the estimated parameter vector
        % theta_est and its Hessian matrix
        [theta_est, ~, ~, ~, ~, ~, hessian] = fmincon(objective_function, ...
            p_initial_vec, [], [], [], [], lb, ub, [], options);
    else
        % Perform optimization to find the estimated parameter vector theta_est
        theta_est = fmincon(objective_function, p_initial_vec, [], [], ...
            [], [], lb, ub, [], options);        
    end
end