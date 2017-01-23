function out = normest(f)
%NORMEST   Estimate the Inf-norm of a CLASSICFUN
%   NORMEST(F) is an estimate of the Inf-norm of the CLASSICFUN F.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = normest(f.onefun);

end
