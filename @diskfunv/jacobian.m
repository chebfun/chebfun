function J = jacobian( F )
%JACOBIAN   Jacobian determinant of a DISKFUNV.
%   J = JACOBIAN(F) computes the determinant of the Jacobian matrix associated
%   with the vector-valued DISKFUNV F. 
%
%   Note we return the determinant of the Jacobian matrix and not the Jacobian
%   matrix itself.
%
% See also DISKFUN/GRADIENT. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( F ) )  
    J = diskfun; 
    return 
end

% Compute derivatives of F: 
Fx = diff( F, 1 ); 
Fy = diff( F, 2 ); 

% The Jacobian is the determinant of the Jacobian matrix: 
Fxc = Fx.components; 
Fyc = Fy.components;
J = Fxc{1} .* Fyc{2} - Fyc{1} .* Fxc{2}; 

end