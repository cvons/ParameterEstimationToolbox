function [g, C] = measurement_fun_3d(x)
    g = x(3);

    if nargout == 2
        C = [0 0 1];
    end
end