function f = abs(f, varargin)
%ABS   Absolute value of a CLASSICFUN object.
%   ABS(F) returns the absolute value of F, where F is a CLASSICFUN object with no
%   roots in F.domain. If ~isempty(roots(F)), then ABS(F) will return garbage
%   with no warning. F may be complex.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Take the absolute value of the ONEFUN:
f.onefun = abs(f.onefun, varargin{:});

end
