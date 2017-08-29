function out = normest(f)
%NORMEST   Estimate the norm of a TRIGTECH.
%   NORMEST(F) estimates the inf norm of a TRIGTECH F simply by looking at the
%   max norm of its values. Use NORMEST when an approximate norm is acceptable.
%
% See also NORM.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Since vscale(f) = max(abs(f.values)) calling max(vscale(f)) is sufficient:
out = max(vscale(f));
            
end
