function J = jacobian(F)
%JACOBIAN   Determinant of the Jacobian of a CHEBFUN3V object.
%   J = JACOBIAN(F) computes the determinant of the Jacobian matrix 
%   associated with the vector-valued CHEBFUN3V object F. F must have three
%   components.
%
%   Note that we return the determinant of the Jacobian matrix and not the 
%   Jacobian matrix itself.
%
% See also CHEBFUN3/GRADIENT. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

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
Fx = diff(F, 1, 1);
Fy = diff(F, 1, 2);
Fz = diff(F, 1, 3);

% Jacobian: 
Fxc = Fx.components; 
Fyc = Fy.components;
Fzc = Fz.components;

J = Fxc{1} .* (Fyc{2} .* Fzc{3} - Fzc{2} .* Fyc{3}) ...
  - Fyc{1} .* (Fxc{2} .* Fzc{3} - Fzc{2} .* Fxc{3}) ...
  + Fzc{1} .* (Fxc{2} .* Fyc{3} - Fyc{2} .* Fxc{3});

end