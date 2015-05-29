function H = cross( F, G )
%CROSS   Vector cross product.
%   CROSS(F, G) returns the cross product of the CHEBFUN2V objects F and G. If F
%   and G both have two components, then it returns the CHEBFUN2 representing
%       CROSS(F,G) = F(1) * G(2) - F(2) * G(1)
%   where F = (F(1); F(2)) and G = (G(1); G(2)). If F and G have three
%   components then it returns the CHEBFUN2V representing the 3D cross 
%   product.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

% Empty check: 
if ( isempty( F ) || isempty( G ) )
    H = chebfun2v;
    return
end

% Get number of components: 
Fc = F.components; 
Gc = G.components; 

% Do curl: 
if ( F.nComponents == 2 && G.nComponents == 2 )  % 2D curl 
    H = Fc{1} .* Gc{2} - Fc{2} .* Gc{1}; 
elseif ( F.nComponents == 3 && G.nComponents == 3 ) % 3D curl
    H = [ Fc{2} .* Gc{3} - Fc{3} .* Gc{2} ; ...
          Fc{3} .* Gc{1} - Fc{1} .* Gc{3} ; ...
          Fc{1} .* Gc{2} - Fc{2} .* Gc{1} ];
else
    error('CHEBFUN:CHEBFUN2V:cross:components', ...
        'CHEBFUN2V objects must be both 2- or 3-vectors.');
end

end
