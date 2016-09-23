function f = dot(F, G)
%DOT   Vector dot product.
%   DOT(F, G) returns the dot product of the SPHEREFUN objects F and G. DOT(F,
%   G) is the same as F'*G.
% 
% See also SPHEREFUNV/CROSS. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

if ( isempty(F) || isempty(G) ) 
    f = spherefun();
    return
end

Fc = F.components; 
Gc = G.components;

f = Fc{1}.*Gc{1} + Fc{2}.*Gc{2} + Fc{3}.*Gc{3};

end
