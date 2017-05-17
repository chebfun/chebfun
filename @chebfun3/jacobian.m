function J = jacobian(f, g, h)
%JACOBIAN   Jacobian determinant of three CHEBFUN3 objects.
%   J = JACOBIAN(F,G,H) returns the determinant of the Jacobian matrix.
%
%   Note we return the determinant of the Jacobian matrix and not the 
%   Jacobian matrix itself.
%
% See also CHEBFUN3V/JACOBIAN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(f) || isempty(g) || isempty(h))  
    J = [];
    return
end

% Call CHEBFUN3V/JACOBIAN():
J = jacobian(chebfun3v({f, g, h}));

end