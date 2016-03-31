function G = gradient( f ) 
%GRADIENT   Numerical surface gradient of a SPHEREFUN. 
%   G = GRADIENT(F) returns the numerical surface gradient of the
%   SPHEREFUN F as a SPHEREFUNV G.
%
% See also SPHEREFUNV/DIV, SPHEREFUNV/CURL, SPHEREFUNV/VORTICITY.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% Empty check.
if ( isempty( f ) )
    G = spherefunv;
    return
end

fx = diff(f, 1);   % diff in x-variable
fy = diff(f, 2);   % diff in y-variable 
fz = diff(f, 3);   % diff in z-variable

G = spherefunv(fx, fy, fz);

end

