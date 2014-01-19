function H = cross( F, G )
%CROSS  Vector cross product.
%
% f = cross(F,G) returns the cross product of the chebfun2v objects
% F and G. If F and G both have two components, then it returns the chebfun2 
% representing
% 
%   CROSS(F,G) = F(1) * G(2) - F(2) * G(1)
% 
% where F = (F(1); F(2)) and G = (G(1); G(2)).  If F and G have three
% components then it returns the chebfun2v representing the 3D cross
% product. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information. 

% Empty check: 
if ( isempty( F ) || isempty( G ) )
    H = chebfun2v;
    return
end

Fc = F.components; 
Gc = G.components; 
if ( F.nComponents == 2 && G.nComponents == 2 )  % 2D curl 
    
    H = Fc(1) .* Gc(2) - Fc(2) .* Gc(1); 
    
elseif ( F.nComponents == 3 && G.nComponents == 3 ) 
    
    H = [ Fc(2).* Gcomp(3) - Fc(3) .* Gc(2);...
      Fc(3) .* Gc(1) - Fc(1) .* Gc(3);...
      Fc(1) .* Gc(2) - Fc(2) .* Gc(1) ];
  
else
    error('CHEBFUN2V:CROSS','Chebfun2v objects must be both 2- or 3-vectors.');
end

end