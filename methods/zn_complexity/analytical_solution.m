% author: mferhata ~02Jan21
function [field, adaptive_rhos, distances_to_shape_center] = analytical_solution (S, s, dt)
    if  ~exist('s', 'var')
        s   = 1;
    end
    if ndims (S) == 3
        S   = padarray(S, [1 1 1]);
    else
        S   = padarray(S, [1 1]);
    end

    if ~exist ('dt', 'var')
        distances_to_boundaries ...
            = bwdist (~S, 'ch');
    else
        distances_to_boundaries = padarray (dt, [1 1 1]);
    end

    rho ...
        = max (distances_to_boundaries(:)) / s;

    exponents ...
        = 1 - distances_to_boundaries / rho;

    field ...
        = rho^2 .* ...
            ( ...
                1 ...
                - ...
                (exp(1) / (exp(2 * s) + 1)) * ...
                (exp (exponents)  + exp (-exponents)) ...
            );

    field(~S)   = 0;
    field       = single (field);
    if ndims(S) == 3
        field       = field(2:end-1, 2:end-1, 2:end-1);
    else
        field       = field(2:end-1, 2:end-1);
    end
end
