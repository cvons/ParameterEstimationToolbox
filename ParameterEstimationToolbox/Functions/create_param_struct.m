function p = create_param_struct(theta, theta_fields, p_fixed)
% create_param_struct Creates a parameter struct from provided inputs.
%
% This function creates a parameter struct using the provided parameter vector
% (theta), field names for the parameters (theta_fields) and a struct of
% fixed parameters (p_fixed). It returns the parameter struct (p).
%
% Inputs:
%   theta: Parameter vector
%   theta_fields: Field names of the parameters
%   p_fixed: Struct of fixed parameters
%
% Output:
%   p: Parameter struct

p = struct;

% Populate parameter struct with values from theta
for i = 1:length(theta)
    p.(theta_fields{i}) = theta(i);
end

% Add fixed parameters to the parameter struct
p_fixed_fields = fieldnames(p_fixed);
for i = 1:length(p_fixed_fields)
    p.(p_fixed_fields{i}) = p_fixed.(p_fixed_fields{i});
end
