function F = uminus( F )
%- Unary minus of a chebfun2v
%   -F returns the chebfun2v negated componentwise. 
%
%   UMINUS(F) is called by the syntax -F. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check: 
if (isempty( F ) )
    return
end

% Take uminus of each component:
F.components = cellfun( @uminus, F.components, 'UniformOutput', false );

end