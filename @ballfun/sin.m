function g = sin(f)
%SIN  Sine of a BALLFUN.
%   SIN(F) computes the sine of the BALLFUN F.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.coeffs;
G = ballfun.vals2coeffs(sin(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
