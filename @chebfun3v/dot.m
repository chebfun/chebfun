function f = dot( F, G )
%DOT   Vector dot product.
%   DOT(F, G) returns the dot product of CHEBFUN3V objects F and G. 
%   DOT(F, G) is the same as F'*G.
% 
%   See also CROSS. 

if ( isempty( F ) || isempty( G ) ) 
    f = chebfun3();
    return
end

nF = F.nComponents; 
nG = G.nComponents; 
if ( nG ~= nF ) 
    error('CHEBFUN:CHEBFUN3V:dot:components', ...
        'CHEBFUN3V objects should have the same number of components.');
end

Fc = F.components; 
Gc = G.components; 
for jj = 1:nF 
    Fc{jj} = times(Fc{jj}, Gc{jj}); 
end

f = chebfun3( 0, G.components{1}.domain );
for jj = 1 : nF
    f = f + Fc{jj};
end

end