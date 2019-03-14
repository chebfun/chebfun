function nm = norm(f)
%NORM  Norm of a BALLFUN.
%   For BALLFUN objects:
%    NORM(F) = sqrt(integral of abs(F)^2).
%   If F is near zero, this function might be inaccurate.
%
%   See also BALLFUNV/NORM.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Do something faster for the test

% Empty check:
if ( isempty(f) )
    nm = [];
    return
end

% Get the coefficients of f
[m,n,p] = size(f);
F = coeffs3(f,2*m,2*n,2*p);

% Convert to values and conpute f*conj(f)
Vals = ballfun.coeffs2vals(F);
G = Vals.*conj(Vals);
g = ballfun.vals2coeffs(G);
g = ballfun(g,'coeffs');

% Compute the norm
nm = sqrt(abs(sum3(g)));
end
