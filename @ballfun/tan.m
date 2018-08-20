function g = tan(f)
% TAN Tangent of a BALLFUN function
%   TAN(f) is the tangent of the BALLFUN function f

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.coeffs;
G = ballfun.vals2coeffs(tan(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
