function J = jacobian( F )
%JACOBIAN   Jacobian determinant of a CHEBFUN3V.
%   J = JACOBIAN(F) computes the determinant of the Jacobian matrix associated
%   to the vector-valued CHEBFUN3V F. The CHEBFUN3V must have three components.
%
%   Note we return the determinant of the Jacobian matrix and not the Jacobian
%   matrix itself.
%
% See also CHEBFUN3/GRADIENT. 

% Empty check. 
if ( isempty(F) )  
    J = []; 
    return 
end

if ( F.nComponents ~= 3 )
    error('CHEBFUN:CHEBFUN3V:jacobian:notSquare', ...
        'Jacobian matrix is not square.')
end

% Determinant formula: 
Fx = diff( F, 1, 1 ); 
Fy = diff( F, 1, 2 ); 
Fz = diff( F, 1, 3 ); 

% Jacobian: 
Fxc = Fx.components; 
Fyc = Fy.components;
Fzc = Fz.components;

J = Fxc{1} .* ( Fyc{2} .* Fzc{3} - Fzc{2} .* Fyc{3} ) ...
  - Fyc{1} .* ( Fxc{2} .* Fzc{3} - Fzc{2} .* Fxc{3} ) ...
  + Fzc{1} .* ( Fxc{2} .* Fyc{3} - Fyc{2} .* Fxc{3} );

end
