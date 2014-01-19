function F = uminus( F )
%- Unary minus of a chebfun2v
% 
% -F returns the chebfun2v negated componentwise. 
% UMINUS(F) is called by the syntax -F. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check: 
if (isempty( F ) )
    return;
end

nF = F.nComponents; 
for jj = 1:nF
    F.components{jj} = uminus( F.components{jj} );
end

end