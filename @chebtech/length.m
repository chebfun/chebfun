function len = length(f)
%LENGTH   Length of a CHEBTECH.
%   LENGTH(F) is the number of values at Chebyshev points used to represent F.
%
% See also SIZE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% The length of a CHEBTECH is the length of its vector of values.
len = size(f.coeffs, 1);

end
