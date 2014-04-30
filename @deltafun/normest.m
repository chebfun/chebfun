function out = normest(f)
%NORMEST   Estimate the norm of a DELTAFUN.
%   NORMEST(F) returns NORMEST of the funPart of a DELTAFUN F.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

out = normest(f.funPart);
            
end
