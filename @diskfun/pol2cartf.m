function fdf = pol2cartf(f,th,r)
% Wrapper for evaluating a function defined in Cartesian coordinates as a
% function defined in polar coordinates.

    x = r.*cos(th);
    y = r.*sin(th);
    
fdf = f(x,y);

end