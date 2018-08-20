function g = sqrt(f)
% SQRT Square root of a BALLFUN function
%   SQRT(f) is the square root of the BALLFUN function f

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.coeffs;
G = ballfun.vals2coeffs(sqrt(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
