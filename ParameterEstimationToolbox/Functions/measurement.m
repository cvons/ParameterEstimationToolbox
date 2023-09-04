function y = measurement(g, X, R, seed)
% measurement Adds measurement noise to the output of a function.
%
% Inputs:
%   g: Function handle to evaluate the input data X
%   X: Input data matrix
%   R: Covariance matrix for the measurement noise
%   seed: Seed for the random number generator (optional)
%
% Outputs:
%   y: Output data with added measurement noise
    
    temp = zeros(numel(g(X(:,1))), size(X,2));
    for i = 1:size(X,2)
        temp(:,i) = g(X(:,i));
    end
    [m, n] = size(temp); % Get the dimensions of temp
    
    if nargin == 4
        v = measurement_noise([m, n], R, seed);
    else
        v = measurement_noise([m, n], R);
    end

    y = temp + v;
end
