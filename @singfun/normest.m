function out = normest(f)
%NORMEST   Estimate the norm of a SINGFUN.
%   NORMEST(F) returns the NORMEST of the smooth part of a SINGFUN F. Since the 
%   function value of a SINGFUN is infinite at -1 and 1 in most of the cases due
%   to the pole(s), NORMEST turn the estimated norm of the smooth part of F.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Call NORMEST() of the smooth part:
out = normest(f.smoothPart);
            
end
