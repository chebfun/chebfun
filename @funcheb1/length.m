function len = length(f)
%LENGTH	  Length of a FUNCHEB1.
%   LENGTH(F) is the number of values at Chebyshev points used to represent F.
%
%   See also SIZE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% The length of a FUNCHEB1 is the length of its vector of values.
len = size(f.values, 1);

end
