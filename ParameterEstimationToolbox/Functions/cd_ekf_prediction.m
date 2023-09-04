function [xhat_k_k1, P_k_k1] = cd_ekf_prediction(fun_jac, fun_diff, ...
    xhat_k1_k1, P_k1_k1, u_k1, p, t_k1, t_k)
% cd_ekf_prediction Performs the prediction step of the CD-EKF.
%
% This function performs the prediction step of the continuous-discrete
% extended Kalman filter (CD-EKF) using the provided drift function and its
% Jacobian (fun_jac), current state estimate (xhat_k1_k1), current state
% covariance matrix (P_k1_k1), previous manipulated variable (u_k1),
% parameter struct (p) and time points (t_k1 and t_k). It returns the predicted
% state estimate (xhat_k_k1) and predicted state covariance matrix (P_k_k1).
%
% Inputs:
%   fun_jac: Drift function and Jacobian handle
%   xhat_k1_k1: Current state estimate
%   P_k1_k1: Current state covariance matrix
%   u_k1: Previous manipulated variable
%   p: Parameter struct
%   t_k1: Previous time point
%   t_k: Current time point
%
% Outputs:
%   xhat_k_k1: Predicted state estimate
%   P_k_k1: Predicted state covariance matrix


N = 1; % Number of steps for prediction of xhat_k_k1
dt = t_k - t_k1; % Time step size

% Compute the one-step prediction of xhat_k_k1 using implicit Euler method
[~, xhat_k_k1] = implicit_euler(fun_jac, t_k1, t_k, ...
    N, xhat_k1_k1, u_k1, p);



% Perform N steps of explicit Euler to predict P_k_k1
P_k_k1 = P_k1_k1;
for i = 1:N
    % Compute A_k1 and sigma_k1 at the current time step
    [~, A_k1] = feval(fun_jac, t_k1+(1-i)*dt, xhat_k_k1(i), u_k1, p);
    sigma_k1 = fun_diff(t_k1+(1-i)*dt, xhat_k_k1(i), u_k1, p);
    
    % Discrete difference equation for P_{k|k-1}
    P_k_k1 = P_k_k1 + dt * (A_k1 * P_k_k1 + P_k_k1 * A_k1' + sigma_k1 * sigma_k1');
end
xhat_k_k1 = xhat_k_k1(:, end);

end