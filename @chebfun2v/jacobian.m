function J = jacobian( F )
%JACOBIAN   Jacobian determinant of a CHEBFUN2V.
%   J = JACOBIAN(F) computes the determinant of the Jacobian matrix associated
%   to the vector-valued CHEBFUN2V F. The CHEBFUN2V must have two components.
%
%   Note we return the determinant of the Jacobian matrix and not the Jacobian
%   matrix itself.
%
% See also CHEBFUN2/GRADIENT. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check. 
if ( isempty(F) )  
    J = []; 
    return 
end

if ( F.nComponents == 3 )
    error('CHEBFUN:CHEBFUN2V:jacobian:notSquare', ...
        'Jacobian matrix is not square.')
end

% Determinant formula: 
Fx = diff( F, 1, 2 ); 
Fy = diff( F, 1, 1 ); 

% Jacobian: 
Fxc = Fx.components; 
Fyc = Fy.components;
J = Fxc{1} .* Fyc{2} - Fyc{1} .* Fxc{2}; 

end
