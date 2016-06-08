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


%   surface curl of F is curl(0, 0,F) = (F_y, -F_x)

G = diskfunv(diff(f, 2), -diff(f, 1));

end

