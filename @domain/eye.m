function I = eye(d)
%EYE       Identity operator.
%   This function is deprecated. Use OPERATORBLOCK.EYE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

I = linop( operatorBlock.eye(d) );

end
