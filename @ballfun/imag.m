function g = imag(f)
%IMAG   Complex imaginary part of a BALLFUN.
%   IMAG(F) is the imaginary part of F.
%
% See also REAL.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.coeffs;
% Compute the imaginary part of the values and return the corresponding array of
% coefficients
G = ballfun.vals2coeffs(imag(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
