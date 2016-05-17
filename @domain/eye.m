function I = eye(d)
%EYE       Identity operator.
%   This function is deprecated. Use OPERATORBLOCK.EYE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

I = linop( operatorBlock.eye(double(d)) );

end
