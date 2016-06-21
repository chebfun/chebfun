function H = cross(F, G)
%CROSS   Cross product of CHEBFUN3V objects.
%   CROSS(F, G) returns a CHEBFUN3V representing the 3D cross product of 
%   the CHEBFUN3V objects F and G. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(F) || isempty(G) )
    H = chebfun3v;
    return
end

% Get number of components:
Fc = F.components; 
Gc = G.components; 

if ( F.nComponents == 3 && G.nComponents == 3 )
    H = [ Fc{2} .* Gc{3} - Fc{3} .* Gc{2};
          Fc{3} .* Gc{1} - Fc{1} .* Gc{3};
          Fc{1} .* Gc{2} - Fc{2} .* Gc{1}];
else
    error('CHEBFUN:CHEBFUN3V:cross:components', ...
        'CHEBFUN3V objects must be both 3-vectors.');
end

end