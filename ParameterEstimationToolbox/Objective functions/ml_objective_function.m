function obj_fun = ml_objective_function(fun_jac, fun_diff, theta, y, U, ...
    x0, P0, T, meas_fun_jac, p_fixed, theta_fields)
% ml_objective_function Computes the negative log-likelihood objective function.
%
% This function calculates the negative log-likelihood objective function
% using the provided drift function and its Jacobian (fun_jac), parameter vector (theta),
% measurement data (y), input trajectory (U), initial state (x0), prior covariance
% (P0), sampling interval (T), measurement function handle including its Jacobian (meas_fun_jac),
% fixed parameters (p_fixed) and field names of the parameters (theta_fields).
%
% Inputs:
%   fun_jac: Drift function handle including its Jacobian
%   theta: Parameter vector
%   y: Measurement data
%   U: Input trajectory
%   x0: Initial state
%   P0: Prior covariance
%   T: Sampling interval
%   meas_fun_jac: Measurement function handle including its Jacobian
%   p_fixed: Struct of fixed parameters
%   theta_fields: Field names of the parameters
%
% Output:
%   obj_fun: Negative log-likelihood objective function value

    % Create the parameter struct from the theta vector and p_fixed
    p = create_param_struct(theta, theta_fields, p_fixed);

    % Initialize CD-EKF
    [~, ~, E_k, R_E_k] = cd_ekf(fun_jac, fun_diff, y, U, x0, P0, T, meas_fun_jac, p);

    % Calculate the negative log-likelihood
    N = size(y, 2)-1;
    n_y = size(y, 1);
    obj_fun = 0;

    for k = 1:N
        obj_fun = obj_fun + log(det(R_E_k(:,:,k))) + E_k(:,k)' * inv(R_E_k(:,:,k)) * E_k(:,k);
    end
    obj_fun = 0.5 * obj_fun + (N + 1) / 2 * n_y * log(2 * pi); 
end