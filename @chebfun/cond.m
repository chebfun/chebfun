function out = cond(f)
%COND   Condition number of a CHEBFUN.
%   COND(F) is the 2-norm condition number of F.
%
% See also SVD, RANK.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Compute the singular values:
s = svd(f, 0);

% Compute the condition number:
out = s(1)/s(end);

end
