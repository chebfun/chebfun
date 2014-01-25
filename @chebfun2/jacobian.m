function J = jacobian( f, g )
%JACOBIAN Jacobian determinant of two chebfun2.
%
% J = JACOBIAN(F,G) returns the Jacobian determinant of the Jacobian
% matrix. 
%
% Note we return the determinant of the Jacobian matrix and not the
% Jacobian matrix itself. 
%
% See also CHEBFUN2V/JACOBIAN. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check: 
if ( isempty( f ) || isempty( g ) )  
    J = []; 
    return
end

% Call chebfun2v/jacobian:
J = jacobian( chebfun2v( {f, g} ) );

end