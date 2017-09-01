function out = normest(f)
%NORMEST   Estimate the Inf-norm of a FUN
%   NORMEST(F) is an estimate of the Inf-norm of the FUN F.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Call normest of the ONEFUN of F
out = normest(f.onefun);

end
