function g = imag(f)
% IMAG Imaginary part of a BALLFUN function
%   IMAG(f) is the imaginary part of the BALLFUN function f

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.coeffs;
% Compute the imaginary part of the values and return the corresponding array of
% coefficients
G = ballfun.vals2coeffs(imag(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
