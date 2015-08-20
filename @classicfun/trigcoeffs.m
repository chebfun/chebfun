function out = trigcoeffs(f, N)
%TRIGCOEFFS   Trigonometric Fourier coefficients of a CLASSICFUN.
%   C = TRIGCOEFFS(F) returns the trigonometric Fourier coefficients of F
%   using complex-exponential form.  Specifically, for N = length(F)
%   If N is odd
%       F(x) = C(1)*z^(N-1)/2 + C(2)*z^((N-1)/2-1) + ... + C((N+1)/2) + ... 
%                + C(N)*z^(-(N-1)/2)
%   If N is even
%       F(x) = C(1)*z^(N/2-1) + C(2)*z^(N/2-2) + ... + C(N/2) + ...
%                + C(N-1)*z^(-N/2-1) + 1/2*C(N)*(z^(N/2) + z^(-N/2))
%   where z = exp(1i*pi*x).
%
%   A = TRIGCOEFFS(F, N) truncates or pads the vector C so that N coefficients
%   of F are returned.
%
%   If F is array-valued with M columns, then C is an MxN matrix.
%
% See also LEGCOEFFS, CHEBCOEFFS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    N = length(f);
end

% Call TRIGCOEFFS() of the .ONEFUN:
out = trigcoeffs(f.onefun, N);

end
