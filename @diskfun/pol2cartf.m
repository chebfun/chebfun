% Wrapper for evaluating a function defined in Cartesian cooridnates as a
% function defined in polar coordinates.
function fdf = pol2cartf(f,th,r)

    x = r.*cos(th);
    y = r.*sin(th);
    

fdf = f(x,y);

end