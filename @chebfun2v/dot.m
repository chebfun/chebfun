function f = dot( F, G )
%DOT   Vector dot product.
%   DOT(F, G) returns the dot product of the CHEBFUN2V objects F and G. DOT(F,
%   G) is the same as F'*G.
% 
% See also CROSS. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

if ( isempty( F ) || isempty( G ) ) 
    f = chebfun2();
    return
end

nF = F.nComponents; 
nG = G.nComponents; 
if ( nG ~= nF ) 
    error('CHEBFUN:CHEBFUN2V:dot:components', ...
        'CHEBFUN2V object should have the same number of components.');
end

Fc = F.components; 
Gc = G.components; 
for jj = 1:nF 
    Fc{jj} = times(Fc{jj}, Gc{jj}); 
end

f = chebfun2( 0, G.components{1}.domain );
for jj = 1 : nF
    f = f + Fc{jj};
end

end
