function f = dot( F, G )
%DOT   Vector dot product.
%   DOT(F, G) returns the dot product of the DISKFUNV objects F and G. 
%   DOT(F, G) is the same as F'*G.
% 
% See also CROSS. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

% Empty check: 
if ( isempty( F ) || isempty( G ) ) 
    f = diskfun();
    return
end

% Extract components: 
Fc = F.components; 
Gc = G.components; 

% Calculate dot-product:
f = times(Fc{1}, Gc{1})+times(Fc{2}, Gc{2}); 

end