function g = exp(f)
% EXP Absolue value of a BALLFUN function
%   EXP(f) is the absolute value of the BALLFUN function f

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.coeffs;
G = ballfun.vals2coeffs(exp(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
