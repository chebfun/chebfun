function J = jacobian( F )
%JACOBIAN   Jacobian determinant of a DISKFUNV.
%   J = JACOBIAN(F) computes the determinant of the Jacobian matrix associated
%   to the vector-valued DISKFUNV F. The DISKFUNV must have two components.
%
%   Note we return the determinant of the Jacobian matrix and not the Jacobian
%   matrix itself.
%
% See also DISKFUN/GRADIENT. 

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check. 
if ( isempty(F) )  
    J = []; 
    return 
end

%if ( F.nComponents == 3 )
 %   error('DISKFUN:DISKFUNV:jacobian:notSquare', ...
  %      'Jacobian matrix is not square.')
%end

% Determinant formula: 
Fx = diff( F,  1 ); 
Fy = diff( F, 2 ); 

% Jacobian: 
Fxc = Fx.components; 
Fyc = Fy.components;
J = Fxc{1} .* Fyc{2} - Fyc{1} .* Fxc{2}; 

end
