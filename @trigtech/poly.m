function out = poly(f)
%POLY   Polynomial coefficients of a TRIGTECH.
%   C = POLY(F) returns the polynomial coefficients of F so that:
%   If N is odd
%       F(x) = C(1)*z^(N-1)/2 + C(2)*z^((N-1)/2-1) + ... + C(N)*z^(-(N-1)/2)
%   If N is even
%       F(x) = C(1)*z^(N/2-1) + C(2)*z^(N/2-2) + ... + C(N-1)*z^(-N/2-1) +
%                  1/2*C(N)*(z^(N/2) + z^(-N/2))
%   where z = exp(1i*pi*x) and -1 <= x <= 1.
%   
%   Note that unlike the MATLAB POLY command, TRIGTECH/POLY can operate on
%   array-valued TRIGTECH objects, and hence produce a matrix output. In such
%   instances, the rows of C correspond to the columns of F = [F1, F2, ...].
%   That is, in the case N is odd
%        F1(x) = C(1,1)*z^(N-1)/2 + C(1,2)*z^((N-1)/2-1) + ... + C(1,N)*z^(-(N-1)/2)
%        F2(x) = C(2,1)*z^(N-1)/2 + C(2,2)*z^((N-1)/2-1) + ... + C(2,N)*z^(-(N-1)/2)
%   This strange behaviour is a result of MATLAB's decision to return a row
%   vector from the POLY command, even for column vector input.
%
% See also CHEBCOEFFS, TRIGCOEFFS, LEGCOEFFS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    out = [];
    return
end

% Simply return the transpose of the coefficients:
out = f.coeffs.';
    
end
