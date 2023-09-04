% Driver Script for Model Simulation and Measurement
% This script demonstrates how to set up and run system simulation and measurement.

% Get the current working directory
currentDir = pwd;

% Add the current directory and all its subfolders to the MATLAB path
addpath(genpath(currentDir));

% Define function handles for the model and measurement functions.
% The function handles are found in \ParameterEstimation\Toolbox\Models
fun_jac_drift = @cstr_1d_fun_jac; % Drift function including Jacobian
fun_diff = @cstr_1d_diffusion;    % Diffusion function
meas_fun = @measurement_fun_1d;  % Measurement function

% Define model parameters
dH = -560;        % Enthalpy change
rho = 1;          % Density
cp = 4.186;       % Heat capacity
sigma_v = 0.15;   % Standard deviation of measurement noise

% Initialize system parameters
p.V = 0.105;            % Reactor volume
p.CAin = 1.6/2;         % Inlet concentration of A
p.CBin = 2.4/2;         % Inlet concentration of B
p.Tin = 273.65;         % Inlet temperature in Kelvin
p.logk0 = 24.6;         % Logarithm of the Arrhenius constant
p.EaR = 8500;           % Activation energy divided by gas constant
p.beta = -dH/(rho*cp);  % Temperature coefficient
p.sigma_T = 5;          % Standard deviation of temperature noise

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

% Run system simulation and measurement
[Y, X, T] = simulate_model(fun_jac_drift, fun_diff, meas_fun, U, [t0, tf], ...
    N, x0, p, seed);

%% Plotting the data
% Create a figure for plotting
fig = figure;
hold on

% Plot the temperature trajectories
plot(T/60, X-273.15, 'LineWidth', 2)

% Add labels to the axes
ylabel('T [\circC]')
xlabel('time [min]')

% Add a legend to the plot
legend('Sim. temp.', 'Location', 'best')