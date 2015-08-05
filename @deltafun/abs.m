function f = abs(f, varargin)
%ABS   Absolute value of a DELTAFUN object.
%   ABS(F) returns the absolute value of F, where F is a DELTAFUN object such
%   that the fun part of F has no roots in the domain of F. If 
%   ~isempty(roots(F)), then ABS(F) will return garbage with no warning.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Take the absolute values:
f.funPart = abs(f.funPart, varargin{:});
f.deltaMag = abs(f.deltaMag);

end
