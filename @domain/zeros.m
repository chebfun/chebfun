function Z = zeros(d)
%ZEROS     Zero operator.
%   This function is deprecated. Use OPERATORBLOCK.ZEROS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

Z = linop( operatorBlock.zeros(double(d)) );

end
