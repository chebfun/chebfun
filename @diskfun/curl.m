function G = curl( f ) 
%CURL   Numerical surface curl of a scalar DISKFUN. 
%   G = CURL(F) returns a DISKFUNV G representing the numerical surface
%   curl of the scalar DISKFUN F. 

% See also GRADIENT.

% Empty check.
if isempty( f )
    G = diskfunv;
    return;
end

% On the disk in Cartesian cooridnates the curl of a scalar is just
% cross(n,gradient(f)), where n is the normal vector (which is just
% [x,y,z]^T).
%   surface curl of F is curl(0, 0,F) = (F_y, -F_x)

G = diskfunv(diff(f, 2), -diff(f, 1));

end

