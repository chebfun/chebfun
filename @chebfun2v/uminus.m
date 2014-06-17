function F = uminus( F )
%- Unary minus of a CHEBFUN2V
%   -F returns the CHEBFUN2V negated componentwise. 
%
%   UMINUS(F) is called by the syntax -F. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if (isempty( F ) )
    return
end

% Take uminus of each component:
F.components = cellfun( @uminus, F.components, 'UniformOutput', false );

end
