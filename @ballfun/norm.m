function nm = norm(f, varargin)
%NORM  Frobenius norm of a BALLFUN.
% 
% For BALLFUN objects:
%    NORM(F) = sqrt(integral of abs(F)^2).

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[m, n, p] = size( f );
if ( m == 1 && n == 1 && p == 1 ) 
    nm = sqrt( 4/3*pi ) * f.coeffs;
else
    nm = abs(sqrt(sum3(f.*conj(f))));
end
end
