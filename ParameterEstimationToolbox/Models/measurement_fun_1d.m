function [g, C] = measurement_fun_1d(x)
    g = x;

    if nargout == 2
        C = 1;
    end
end