function I = integral3(f)
%INTEGRAL3   Triple integral of a BALLFUN over its domain.
%   I = INTEGRAL3(F) returns the double definite integral of a BALLFUN.
%
% See also SUM3.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

I = sum3(f);
end