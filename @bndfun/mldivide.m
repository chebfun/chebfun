function X = mldivide(A, B)
%\   Left matrix divide for BNDFUN objects.
%
%   A\B returns the least-squares solution (with respect to the continuous L^2
%   norm) of A*X = B where A and B are BNDFUN objects.
%
% See also QR, MRDIVIDE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Compute X. Note that since the inner product between functions on [a,b] is
% just a scaling of the inner product between functions on [-1,1], the
% least-squares solution will have the same coefficients as we obtain from
% calling \ on the ONEFUN objects of A and B.
X = A.onefun\B.onefun;

end
