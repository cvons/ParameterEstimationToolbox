function ml_term = map_objective_function(fun_jac, fun_diff, theta, Y, ...
    U, x0, P0, T, meas_fun_jac, p_fixed, p_prior, P_theta_0, theta_fields)
% map_objective_function Computes the MAP (maximum a posteriori) objective function.
%
% This function calculates the MAP (maximum a posteriori) objective function
% using the provided drift function and its Jacobian (fun_jac), parameter vector (theta),
% measurement data (Y), input trajectory (U), initial state (x0), prior covariance
% (P0), time points (T), measurement function handle including its Jacobian (meas_fun_jac),
% fixed parameters (p_fixed), prior parameters (p_prior), prior covariance matrix (P_theta_0)
% and field names of the parameters (theta_fields).
%
% Inputs:
%   fun_jac: Drift function handle including its Jacobian
%   theta: Parameter vector
%   Y: Measurement data
%   U: Input trajectory
%   x0: Initial state
%   P0: Prior covariance
%   T: Time points
%   meas_fun_jac: Measurement function handle including its Jacobian
%   p_fixed: Struct of fixed parameters
%   p_prior: Prior parameter struct
%   P_theta_0: Prior covariance matrix for MAP estimation
%   theta_fields: Field names of the parameters
%
% Output:
%   obj_fun: MAP (maximum a posteriori) objective function value

    % Create the parameter struct from the theta vector and p_fixed
    p = create_param_struct(theta, theta_fields, p_fixed);

    % Initialize CD-EKF
    [~, ~, E_k, R_E_k] = cd_ekf(fun_jac, fun_diff, Y, U, x0, P0, T, meas_fun_jac, p);

    % Calculate the negative log-likelihood (ML objective)
    N = size(Y, 2) - 1;
    n_y = size(Y, 1);
    ml_term = 0;

    for k = 1:N
        ml_term = ml_term + log(det(R_E_k(:,:,k))) + E_k(:,k)' * inv(R_E_k(:,:,k)) * E_k(:,k);
    end
    ml_term = 0.5 * ml_term + (N + 1) / 2 * n_y * log(2 * pi); 

    % Convert the p_prior struct to an array
    p_prior_fields = fieldnames(p_prior);
    theta_prior_vec = zeros(numel(p_prior_fields), 1);
    for i = 1:numel(p_prior_fields)
        theta_prior_vec(i) = p_prior.(p_prior_fields{i});
    end

    % Calculate the MAP regularization term
    theta_diff = theta - theta_prior_vec;
    n_theta = length(theta);
    map_term = 0.5 * (theta_diff' * inv(P_theta_0) * theta_diff + ...
        log(det(P_theta_0)) + n_theta/2 * log(2 * pi));

    % Combine the ML and MAP terms to form the MAP objective
    ml_term = ml_term + map_term;
end