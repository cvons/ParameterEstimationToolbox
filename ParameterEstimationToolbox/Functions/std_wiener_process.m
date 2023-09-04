function [W, Tw, dW] = std_wiener_process(time_interval, N, nW, seed)
% STD_WIENER_PROCESS Generates realizations of a standard Wiener process.
%
% Syntax: [W, Tw, dW] = std_wiener_process(time_interval, N, nW, seed)
% W   : Standard Wiener process in time_interval
% Tw  : Time points
% dW  : White noise used to generate the Wiener process
%
% time_interval: [t0, tf] Time interval for simulation
% N            : Number of intervals
% nW           : Dimension of W(k)
% seed         : Seed for the random number generator (optional)

% Check if seed argument is provided, otherwise use the default random seed
if nargin == 5 && ~isempty(seed)
    rng(seed);
end

t0 = time_interval(1);
tf = time_interval(2);
dt = (tf - t0) / N; % Time step
dW = sqrt(dt) * randn(nW, N); % Generate white noise
W = [zeros(nW, 1), cumsum(dW, 2)]; % Calculate Wiener process
Tw = t0:dt:tf; % Time points
end
