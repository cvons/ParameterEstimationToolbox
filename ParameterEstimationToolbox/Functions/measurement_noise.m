function v = measurement_noise(dim, R, seed)
% measurement_noise Generates measurement noise for a given covariance matrix.
%
% Inputs:
%   dim: Dimensions of the output noise matrix [rows, columns]
%   R: Covariance matrix for the noise
%   seed: Seed for the random number generator (optional)
%
% Outputs:
%   v: Vector of generated measurement noise

    if nargin == 3 && ~isempty(seed)
        rng(seed);
    end
    
    % Generate noise for each column of the matrix
    v = sqrt(diag(R)) .* randn(dim);
end
