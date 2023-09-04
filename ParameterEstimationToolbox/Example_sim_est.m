% Driver Script for Parameter Estimation
% This script demonstrates how to set up and run SDE simulation and 
% parameter estimation.

% Get the current working directory
currentDir = pwd;

% Add the current directory and all its subfolders to the MATLAB path
addpath(genpath(currentDir));

% Define function handles for the model and measurement functions.
% The function handles are found in \ParameterEstimation\Toolbox\Models
fun_jac_drift = @cstr_1d_fun_jac; % Drift function including Jacobian
fun_diff = @cstr_1d_diffusion;    % Diffusion function
meas_fun_jac = @measurement_fun_1d; % Measurement function and its Jacobian

% Define model parameters
dH = -560;        % Enthalpy change
rho = 1;          % Density
cp = 4.186;       % Heat capacity
sigma_v = 0.15;   % Standard deviation of measurement noise

% Initialize system parameters
p.V = 0.105;       % Reactor volume
p.CAin = 1.6/2;    % Inlet concentration of A
p.CBin = 2.4/2;    % Inlet concentration of B
p.Tin = 273.65;    % Inlet temperature in Kelvin
p.logk0 = 24.6;    % Logarithm of the Arrhenius constant
p.EaR = 8500;      % Activation energy divided by gas constant
p.beta = -dH/(rho*cp); % Temperature coefficient
p.sigma_T = 5;     % Standard deviation of temperature noise

% Covariance matrix for measurement noise
% Note: The field named 'R' with the covariance matrix is a requirement for
% the toolbox. While the parameter struct can have any name, it must 
% contain a field named 'R'.
p.R = sigma_v * sigma_v';

% Random seed for reproducibility
seed = 1401;

% Define simulation time
t0 = 0;
tf = 2100;
% Define number of simulation intervals
N = 35000;

% Initial condition for the system
x0 = 273.65;

% Define input trajectory
N_per_min = N/35;
U = [repmat(0.7/60,1,3*N_per_min), repmat(0.6/60,1,2*N_per_min), ...
    repmat(0.5/60,1,2*N_per_min), repmat(0.4/60,1,2*N_per_min), ...
    repmat(0.3/60,1,3*N_per_min), repmat(0.2/60,1,4*N_per_min), ... 
    repmat(0.3/60,1,2*N_per_min), repmat(0.4/60,1,2*N_per_min), ...
    repmat(0.5/60,1,2*N_per_min), repmat(0.6/60,1,2*N_per_min), ...
    repmat(0.7/60,1,4*N_per_min), repmat(0.2/60,1,4*N_per_min), ...
    repmat(0.7/60,1,3*N_per_min)];

% Initial guess for the state vector and covariance matrix
x0_est = x0;
P0_est = 0.0001;

% Initial guesses for parameters
p_initial_guess.beta = 130;
p_initial_guess.logk0 = 24;
p_initial_guess.EaR = 8400;
p_initial_guess.sigma_T = 4;
p_initial_guess.R = 0.1;

% Define parameter bounds
lb = [120,23,8000,1,0];
ub = [140,25.5,9000,10,3];

% Type of objective function for parameter estimation ('ls' for least 
% squares, 'ml' for maximum likelihood, 'map' for maximum a posteriori)
obj_fun_type = 'ml';

% Run parameter estimation
[p_est, Y, X, T, hessian] = param_est_sim(fun_jac_drift, fun_diff, ...
    meas_fun_jac, U, [t0, tf], N, x0, p, p_initial_guess, obj_fun_type, ...
    x0_est, P0_est, lb, ub, seed);


%% Plotting the data
% Update the parameters with the estimated values for simulation
p.beta = p_est(1);
p.logk0 = p_est(2);
p.EaR = p_est(3);
p.sigma_T = p_est(4);
p.R = p_est(5);

% Optionally, set a new random seed for reproducible simulation.
% Note: Choosing a new seed ensures different Wiener process values are
% used for simulation with nominal and estimated parameters. If no new seed
% is selected, the previous seed remains, leading to identical Wiener
% process values for both simulations.
seed = 2203;

% Simulate the system with the estimated parameters
[~, X_est, ~] = simulate_model(fun_jac_drift, fun_diff, ...
    meas_fun_jac, U, [t0, tf], N, x0, p, seed);

% Create a figure for plotting
fig = figure;
hold on

% Plot the original and estimated temperature trajectories
plot(T/60, X-273.15, T/60, X_est-273.15, '--', 'LineWidth', 2)

% Add labels to the axes
ylabel('T [\circC]')
xlabel('time [min]')

% Add a legend to the plot
legend('Nom. param.', 'Est. param.', 'Location', 'best')