function [xhat_k_k, P_k_k, e_k, R_e_k] = cd_ekf_filtering(xhat_k_k1, P_k_k1, y_k, R_k, meas_fun_jac)
% cd_ekf_filtering Performs the filtering step of the CD-EKF.
%
% This function performs the filtering step of the continuous-discrete
% extended Kalman filter (CD-EKF) using the provided current state estimate
% (xhat_k_k1), current state covariance matrix (P_k_k1), measurement (y_k),
% measurement noise covariance matrix (R_k) and measurement function
% with its Jacobian (meas_fun_jac). It returns the updated state estimate (xhat_k_k),
% updated state covariance matrix (P_k_k), innovation (e_k) and innovation
% covariance matrix (R_e_k).
%
% Inputs:
%   xhat_k_k1: Current state estimate
%   P_k_k1: Current state covariance matrix
%   y_k: Measurement at current step
%   R_k: Measurement noise covariance matrix
%   meas_fun_jac: Measurement function and Jacobian handle
%
% Outputs:
%   xhat_k_k: Updated state estimate
%   P_k_k: Updated state covariance matrix
%   e_k: Innovation
%   R_e_k: Innovation covariance matrix

% Compute the measurement estimate and its Jacobian
[yhat_k_k1, C_k] = meas_fun_jac(xhat_k_k1);

% Calculate the innovation and innovation covariance
e_k = y_k - yhat_k_k1;
R_e_k = R_k + C_k * P_k_k1 * C_k';

% Compute the Kalman gain
K_k = P_k_k1 * C_k' * inv(R_e_k);

% Update the state estimate and its covariance
xhat_k_k = xhat_k_k1 + K_k * e_k;
P_k_k = (eye(size(P_k_k1)) - K_k * C_k) * P_k_k1 * (eye(size(P_k_k1)) - K_k * C_k)' + K_k * R_k * K_k';
