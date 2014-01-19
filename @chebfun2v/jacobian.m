function J = jacobian( F )
%JACOBIAN Jacobian determinant of a chebfun2v.
%
% J = JACOBIAN(F) computes the determinant of the Jacobian matrix
% associated to the vector-valued chebfun2v F. The chebfun2v must have two
% components. 
%
% Note we return the determinant of the Jacobian matrix and not the
% Jacobian matrix itself. 
%
% See also CHEBFUN2/GRADIENT. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isempty(F) )  % empty check. 
    J = []; 
    return 
end

if ( F.nComponents == 3 )
    error('CHEBFUN2V:JACOBIAN','Jacobian matrix is not square.')
end

% Determinant formula: 
Fx = diff( F, 1, 2); 
Fy = diff( F, 1, 1); 

% Jacobian: 
J = Fx.components{1}.*Fy.components{2} - Fy.components{1}.*Fx.components{2}; 

end
