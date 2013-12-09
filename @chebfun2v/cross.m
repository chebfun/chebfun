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


Fcomp = F.components; 
Gcomp = G.components; 
if ( F.nComponents == 2 && G.nComponents == 2 )  % no third component then return the 2D curl. 
    H = Fcomp(1) .* Gcomp(2) - Fcomp(2) .* Gcomp(1); % chebfun2 object. 
elseif ( F.nComponents == 3 && G.nComponents == 3 ) 
    H = [ Fcomp(2).* Gcomp(3) - Fcomp(3) .* Gcomp(2);...
      Fcomp(3) .* Gcomp(1) - Fcomp(1) .* Gcomp(3);...
      Fcomp(1) .* Gcomp(2) - Fcomp(2) .* Gcomp(1) ];
else
    error('CHEBFUN2V:CROSS','Chebfun2v objects must be both 2- or 3-vectors.');
end

end