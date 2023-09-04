function ns = struct_subtract(s1, s2)
% struct_subtract Creates a new struct by subtracting fields from another struct.
%
% This function creates a new struct by copying all the fields from s1 that
% are not present in s2. It returns the new struct with the selected fields.
%
% Inputs:
%   s1: Source struct from which fields will be copied
%   s2: Struct containing fields to exclude from the new struct
%
% Output:
%   ns: New struct containing fields from s1 that are not in s2

s1_fields = fieldnames(s1);
s2_fields = fieldnames(s2);

% Initialize the new struct
ns = struct();

% Iterate through the fields of s1
for i = 1:numel(s1_fields)
    field_name = s1_fields{i};
    
    % Check if the field is not present in s2
    if ~any(strcmp(field_name, s2_fields))
        % Add the field and its value to ns
        ns.(field_name) = s1.(field_name);
    end
end