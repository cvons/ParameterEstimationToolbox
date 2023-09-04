% Driver Script for Parameter Estimation from saved data. This can either
% be from an earlier simulation or from experiments.
% This script demonstrates how to load measured data, input trajectory, and
% time points for measurements from a saved data file (of type .mat) and 
% perform parameter estimation.

% Load the saved data file
load('SavedData.mat');  % Make sure SavedData.mat contains variables Y, U, and T

% Define function handles for the model and measurement functions.
% The function handles are found in \ParameterEstimation\Toolbox\Models
fun_jac_drift = @cstr_1d_fun_jac; % Drift function including Jacobian
fun_diff = @cstr_1d_diffusion;    % Diffusion function
meas_fun_jac = @measurement_fun_1d; % Measurement function and its Jacobian

% Define system parameters.
% Note: Include all system parameters, including those to be estimated, 
% in the 'p' struct. Initialization values for parameters to be estimated 
% are arbitrary (for the 'p' struct), as they are determined during the 
% estimation process.
p.V = 0.105;       % Reactor volume
p.CAin = 1.6/2;    % Inlet concentration of A
p.CBin = 2.4/2;    % Inlet concentration of B
p.Tin = 273.65;    % Inlet temperature in Kelvin
p.logk0 = 0;       % Logarithm of the Arrhenius constant
p.EaR = 0;         % Activation energy divided by gas constant
p.beta = 0;        % Temperature coefficient
p.sigma_T = 0;     % Standard deviation of temperature noise
p.R = 0;           % Covariance matrix for measurement noise

% Define initial guesses for parameters (p_guess) and bounds (lb and ub)
p_initial_guess.beta = 130;
p_initial_guess.logk0 = 24;
p_initial_guess.EaR = 8400;
p_initial_guess.sigma_T = 4;
p_initial_guess.R = 0.1;

lb = [120, 23, 8000, 1, 0];
ub = [140, 25.5, 9000, 10, 3];

% Define other required inputs
x0_est = 273.65;     % Initial guess for the state vector
P0_est = 0.0001;     % Initial covariance matrix for the state vector
obj_fun_type = 'ml'; % Type of objective function for parameter estimation

% Run parameter estimation using the loaded data
[p_est, hessian] = estimate_parameters(fun_jac_drift, fun_diff, Y, U, x0_est, P0_est, ...
    T, p, p_initial_guess, meas_fun_jac, obj_fun_type, lb, ub);

%% Plotting the data
% Update the parameters with the estimated values for simulation
p.beta = p_est(1);
p.logk0 = p_est(2);
p.EaR = p_est(3);
p.sigma_T = p_est(4);
p.R = p_est(5);

% Optionally, set a new random seed for reproducible simulation.
% Note: The loaded data was simulated using seed 1401.
% Choosing a new seed ensures different Wiener process values are
% used for simulation with nominal and estimated parameters. If 1401 is 
% chosen, it will lead to identical Wiener process values for both simulations.
seed = 2203;

% Simulate the system with the estimated parameters
[Y_est, ~, ~] = simulate_model(fun_jac_drift, fun_diff, ...
    meas_fun_jac, U, [T(1), T(end)], numel(U), x0_est, p, seed);

% Create a figure for plotting
fig = figure;
hold on

% Plot the original and estimated measurement trajectories (Y)
plot(T/60, Y, T/60, Y_est, '--', 'LineWidth', 2)

% Add labels to the axes
ylabel('T [\circC] with measurement noise')
xlabel('time [min]')

% Add a legend to the plot
legend('Nom. param.', 'Est. param.', 'Location', 'best')
