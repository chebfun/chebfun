function H = dot( F, G )
%DOT   Vector dot product.
%   DOT(F, G) returns the dot product of the CHEBFUN2V objects F and G. DOT(F,
%   G) is the same as F'*G.
% 
% See also CROSS. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

if ( isempty( F ) ) 
    H = F;
    return
end

if ( isempty( G ) ) 
    H = G;
    return
end

nF = F.nComponents; 
nG = G.nComponents; 

% Check that F and G have the same number of components:
if ( nG ~= nF ) 
    error('CHEBFUN:CHEBFUN2V:dot:components', ...
        'CHEBFUN2V object should have the same number of components.');
end

Fc = F.components; 
Gc = G.components; 
for jj = 1:nF 
    Fc{jj} = times(Fc{jj}, Gc{jj}); 
end

H = Fc{1};
for jj = 2 : nF
    H = H + Fc{jj};
end

end
