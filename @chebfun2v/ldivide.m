function H = ldivide( F, G )
%.\   Pointwise chebfun2v left divide.
%
% Left componentwise divide for a chebfun2v. 
%
% See also RDIVIDE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( ( isempty(F) ) || ( isempty(G) ) )
   H = chebfun2v;
   return
end

nF = F.ncomponents;
nG = G.ncomponents;
if ( nF ~= nG ) 
    error('CHEBFUN2V:LDIVIDE','Chebfun2v do not have the same number of components.')
end

H = F; 
for j = 1 : F.ncomponents
    H.components{j} = ldivide( F.components{j}, G.components{j} );
end

end