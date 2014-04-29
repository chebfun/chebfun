function X = mrdivide(A, B)
%/   Right matrix divide for SINGFUN objects.
%
%   Since SINGFUN does not support array-valued objects, right matrix divide is 
%   exactly same as right divide.
%
% See also RDIVIDE, LDIVIDE, MLDIVIDE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

X = rdivide(A, B);

end
