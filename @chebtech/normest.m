function out = normest(f)
%NORMEST   Estimate the norm of a CHEBTECH.
%   NORMEST(F) estimates the inf norm of a CHEBTECH F simply by looking at the
%   max norm of its values. Use NORMEST when an approximate norm is acceptable.
%
% See also NORM.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Since vscale(f) = max(abs(f.values)) this is sufficient:
out = max(vscale(f));
            
end
