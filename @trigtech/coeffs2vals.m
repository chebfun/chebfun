function values = coeffs2vals(coeffs)
%COEFFS2VALS   Convert Fourier coefficients to values at N equally spaced
%points between [-1 1), where N is the number of coefficients.
%   V = COEFFS2VALS(C) returns the values of the trignometric polynomial 
%   as follows:
%   If N is odd
%       F(x) = C(1)*z^(-(N-1)/2) + C(2)*z^(-(N-1)/2-1) + ... + C(N)*z^((N-1)/2)
%   If N is even
%       F(x) = C(1)*z^(-N/2) + C(2)*z^(-N/2+1) + ... + C(N)*z^(N/2-1)
%   where z = exp(1i*pi*x) and -1 <= x <= 1. 
%
%   If the input C is an (N+1)xM matrix then V = COEFFS2VALS(C) returns the
%   (N+1)xM matrix of values V such that V(i,j) is the ith value of the
%   trignometric polynomial corresponding to the jth column.
%
% See also VALS2COEFFS, TRIGPTS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Get the length of the input:
n = size(coeffs, 1);

% Trivial case (constant or empty):
if ( n <= 1 )
    values = coeffs; 
    return
end

% The coefficients are for interpolation defined on [-pi,pi), but the FFT
% works for values on [0,2*pi). To fix the coefficients for this we just need to
% assign c_k = (-1)^k c_k, with k=-(N-1)/2:(N-1)/2 for N odd, and 
% k = -N/2:N/2-1 for N even.
if ( mod(n, 2) ) 
    even_odd_fix = (-1).^(-(n-1)/2:(n-1)/2).';
else
    even_odd_fix = (-1).^((-n/2):(n/2-1)).';
end
coeffs = bsxfun(@times, coeffs, even_odd_fix);

% Shift the coefficients properly.
values = ifft(ifftshift(n*coeffs, 1), [], 1);

end
