function n = norm(f, varargin)
%   NORM Norm of a BALLFUN function.
% 
% For BALLFUN objects:
%    NORM(F) = sqrt(integral of abs(F)^2).
%    NORM(F, inf) = global maximum in absolute value. NOT IMPLEMENTED.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

n = abs(sqrt(sum3(abs(f).^2)));
end
