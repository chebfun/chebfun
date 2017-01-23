function M = diag(f)
%DIAG   Multiplication operator.
%   M = DIAG(F) returns a multiplication operator M so that M*u is equivalent to
%   F.*u.
%
%   This method is mainly provided for backwards compatibility.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

M = operatorBlock.mult(f);

end
