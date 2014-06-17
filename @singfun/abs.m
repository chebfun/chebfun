function f = abs(f, varargin)
%ABS   Absolute value of a SINGFUN object.
%   ABS(F) returns the absolute value of F, where F is a SINGFUN object with no
%   roots in [-1 1]. If ~isempty(roots(F)), then ABS(F) will return garbage
%   with no warning. F may be complex.

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

f.smoothPart = abs(f.smoothPart);

end
