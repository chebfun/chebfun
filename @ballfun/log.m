function g = log(f)
% LOG Logarithm of a BALLFUN function
%   LOG(f) is the logarithm of the BALLFUN function f

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Return the logarithm of the ballfun function f
F = f.coeffs;
G = ballfun.vals2coeffs(log(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
