function I = integral(f)
%INTEGRAL   Triple integral of a BALLFUN over its domain.
%   I = INTEGRAL(F) returns the double definite integral of a BALLFUN.
%
% See also SUM3.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

I = sum3(f);
end