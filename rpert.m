%% RPERT
function number = rpert(x_min, x_max, x_mode, lambda)
% generates a random number from a PERT distribution
% inputs: minimum and maximum values, mode
    %lambda = 400; %   default value
    if(x_min > x_max || x_mode > x_max || x_mode < x_min )
        fprintf("invalid parameters");
    end
    x_range = x_max - x_min;
    if (x_range == 0) 
        number = x_min;
    end
    mu = (x_min + x_max + lambda*x_mode)/(lambda + 2);

    % special case if mu == mode
    if( mu == x_mode )
        v = (lambda/2) + 1;
    else
        v = ((mu - x_min)*(2*x_mode - x_min - x_max))/((x_mode - mu)*(x_max - x_min));
    end
    w = ( v * ( x_max - mu ))/( mu - x_min );
    number = random('Beta', v, w)* x_range + x_min;
end
    