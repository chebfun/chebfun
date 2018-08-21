function g = exp(f)
%EXP   Exponential of a BALLFUN.
%   EXP(F) computes the exponential of the BALLFUN F.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.coeffs;
G = ballfun.vals2coeffs(exp(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
