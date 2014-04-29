function S = sum(d)
%SUM       Integration functional.
%   This function is deprecated. Use FUNCTIONALBLOCK.SUM.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

S = linop( functionalBlock.sum(d) );

end
