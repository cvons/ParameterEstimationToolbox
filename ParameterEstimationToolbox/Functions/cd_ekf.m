function [Xhat, P, E_k, R_E_k] = cd_ekf(fun_jac, fun_diff, Y, U, x0, P0, ...
    T, meas_fun_jac, p)
% cd_ekf Implements the continuous-discrete extended Kalman filter (CD-EKF).
%
% This function performs continuous-discrete extended Kalman filtering
% using the provided drift function and its Jacobian (fun_jac), measurement data (Y),
% input trajectory (U), initial state estimate (x0), initial state covariance
% (P0), time points (T), measurement function and its Jacobian (meas_fun_jac),
% and parameter struct (p). It returns the filtered state estimates (Xhat),
% state covariance matrices (P), innovation (E_k) and innovation covariance (R_E_k).
%
% Inputs:
%   fun_jac: Drift function and Jacobian handle
%   Y: Measurement data
%   U: Input trajectory
%   x0: Initial state estimate
%   P0: Initial state covariance
%   T: Time points
%   meas_fun_jac: Measurement function and Jacobian handle
%   p: Parameter struct
%
% Outputs:
%   Xhat: Filtered state estimates
%   P: State covariance matrices
%   E_k: Innovation
%   R_E_k: Innovation covariance

% Get dimensions
nx = length(x0);   % Number of states
ny = size(Y, 1);   % Number of measurements
N = length(T); % Number of time steps

% Initialize arrays for storing results
Xhat = zeros(nx, N);   % State estimates
P = zeros(nx, nx, N);  % State covariance matrices
E_k = zeros(ny, N);        % Innovations
R_E_k = zeros(ny, ny, N);  % Innovation covariances

% Measurement noise covariance
Rk = p.R;

% % Generate initial state as a normally distributed variable with Niid(x0; P0)
% x0 = mvnrnd(x0, P0)';

% Set initial values for xhat_k1_k1 and P_k1_k1
xhat_k_k1 = x0;
P_k_k1 = P0;

y_k = Y(:, 1);

% Initial filtering step
[xhat_k_k, P_k_k, e_k, R_e_k] = cd_ekf_filtering(xhat_k_k1, P_k_k1, ...
    y_k, Rk, meas_fun_jac);

% Update for the next iteration
xhat_k1_k1 = xhat_k_k;
P_k1_k1 = P_k_k;

% Store the results for each iteration
Xhat(:, 1) = xhat_k1_k1;
P(:, :, 1) = P_k1_k1;
E_k(:, 1) = e_k;
R_E_k(:, :, 1) = R_e_k;

% Perform prediction and filtering for each time step
for k = 2:N
    u_k1 = U(k - 1);   % Previous manipulated variable
    y_k = Y(:, k);     % Measurement at current step
    
    % Prediction step
    [xhat_k_k1, P_k_k1] = cd_ekf_prediction(fun_jac, fun_diff, xhat_k1_k1, ...
        P_k1_k1, u_k1, p, T(k - 1), T(k));
    
    % Filtering step
    [xhat_k_k, P_k_k, e_k, R_e_k] = cd_ekf_filtering(xhat_k_k1, P_k_k1, ...
        y_k, Rk, meas_fun_jac);
    
    % Update for the next iteration
    xhat_k1_k1 = xhat_k_k;
    P_k1_k1 = P_k_k;
    
    % Store the results for each iteration
    Xhat(:, k) = xhat_k1_k1;
    P(:, :, k) = P_k1_k1;
    E_k(:, k) = e_k;
    R_E_k(:, :, k) = R_e_k;
end
