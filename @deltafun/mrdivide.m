function X = mrdivide(A, B)
%/   Right matrix divide for DELTAFUN objects.
%
%   DELTAFUN (at the moment) does not support array-valued objects, 
%   right matrix divide is exactly the same as right divide.
%
% See also RDIVIDE, LDIVIDE, MLDIVIDE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

X = rdivide(A, B);

end
