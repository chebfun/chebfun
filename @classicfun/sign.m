function f = sign(f, varargin)
%SIGN   Signum of a CLASSICFUN object.
%   SIGN(F) returns the sign of F, where F is a CLASSICFUN object with no roots in
%   its domain. If F has roots, then SIGN(F) will return garbage with no
%   warning.
%
%   For the nonzero elements of complex F, SIGN(F) = F ./ ABS(F).
%
% See also ABS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Take the absolute value of the ONEFUN:
f.onefun = sign(f.onefun, varargin{:});

end
