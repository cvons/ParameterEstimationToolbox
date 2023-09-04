function obj_fun = ls_objective_function(fun_jac, fun_diff, theta, Y, U, ...
    x0, P0, T, meas_fun_jac, p_fixed, theta_fields)
% ls_objective_function Computes the least squares objective function.
%
% This function calculates the least squares objective function using the
% provided drift function and its Jacobian (fun_jac), parameter vector (theta),
% measurement data (Y), input trajectory (U), initial state (x0), prior covariance
% (P0), time points (T), measurement function and its Jacobian (meas_fun_jac),
% fixed parameters (p_fixed) and field names of the parameters (theta_fields).
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
%   theta_fields: Field names of the parameters
%
% Output:
%   obj_fun: Least squares objective function value

    % Create the parameter struct from the theta vector and p_fixed
    p = create_param_struct(theta, theta_fields, p_fixed);

    % Initialize CD-EKF
    [~, ~, E_k, ~] = cd_ekf(fun_jac, fun_diff, Y, U, x0, P0, T, meas_fun_jac, p);

    % Calculate the least squares objective function
    obj_fun = 0.5 * sum(sum(E_k.^2));
end