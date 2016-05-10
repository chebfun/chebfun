function S = sum(d)
%SUM       Integration functional.
%   This function is deprecated. Use FUNCTIONALBLOCK.SUM.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

S = linop( functionalBlock.sum(double(d)) );

end
