function f = imag(f)
%IMAG   Complex imaginary part of a BALLFUN.
%   IMAG(F) is the imaginary part of F.
%
% See also REAL.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f.coeffs = ballfun.vals2coeffs( imag( ballfun.coeffs2vals( f.coeffs ) ) );
f = simplify(f);
end
