function J = jacobian( f, g )
%JACOBIAN   Jacobian determinant of two CHEBFUN2.
%   J = JACOBIAN(F,G) returns the Jacobian determinant of the Jacobian matrix.
%
%   Note we return the determinant of the Jacobian matrix and not the Jacobian
%   matrix itself.
%
% See also CHEBFUN2V/JACOBIAN. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( f ) || isempty( g ) )  
    J = []; 
    return
end

% Call CHEBFUN2V/JACOBIAN():
J = jacobian( chebfun2v( { f, g } ) );

end