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

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% *Note about symmetry*.  Some of the code below is designed to
% enforce two symmetries whose failure might disturb users:
% COEFFS exactly real ==> VALUES exactly hermitian
% COEFFS exactly imaginary ==> VALUES exactly skew-hermitian
% This is necessary because the MATLAB FFT code does not
% exactly preserve symmetries.

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

% test for symmetry
isHerm = max(abs(imag(coeffs)),[],1) == 0;
isSkew = max(abs(real(coeffs)),[],1) == 0;

% Shift the coefficients properly.
values = ifft(ifftshift(n*coeffs, 1), [], 1);

% correct if symmetric
vals = [values;values(1,:)];
hermvals = (vals+flipud(conj(vals)))/2;
skewvals = (vals-flipud(conj(vals)))/2;
values(:,isHerm) = hermvals(1:end-1,isHerm);
values(:,isSkew) = skewvals(1:end-1,isSkew);

end
