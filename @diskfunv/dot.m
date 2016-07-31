function f = dot( F, G )
%DOT   Vector dot product.
%   DOT(F, G) returns the dot product of the DISKFUNV objects F and G. DOT(F,
%   G) is the same as F'*G.
% 
% See also CROSS. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

if ( isempty( F ) || isempty( G ) ) 
    f = diskfun();
    return
end


Fc = F.components; 
Gc = G.components; 

f = times(Fc{1}, Gc{1})+times(Fc{2}, Gc{2}); 


end
