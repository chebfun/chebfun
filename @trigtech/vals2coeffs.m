function coeffs = vals2coeffs(values)
%VALS2COEFFS   Convert values at N equally spaced points between [-1 1) 
%   to N trigonometric coefficients.
%   C = VALS2COEFFS(V) returns the vector of N coefficients such that: 
%   If N is odd
%       F(x) = C(1)*z^(-(N-1)/2) + C(2)*z^(-(N-1)/2-1) + ... + C(N)*z^((N-1)/2)
%   If N is even
%       F(x) = C(1)*z^(-N/2) + C(2)*z^(-N/2+1) + ... + C(N)*z^(N/2-1)           
%   where z = exp(1i*pi*x) and -1 <= x <= 1. 
%
%   F(x) interpolates the data [V(1) ; ... ; V(N)] at the N equally 
%   spaced points x_k = -1 + 2*k/N, k=0:N-1. 
%
%   If the input V is an (N+1)xM matrix, then C = VALS2COEFFS(V) returns the
%   (N+1)xM matrix of coefficients C such that F_j(x) intterpolates
%   [V(1,j) ; ... ; V(N+1,j)] for j=1:M using the same formula as above for
%   each column.
%
% See also COEFFS2VALS, TRIGPTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% *Note about symmetry*.  Some of the code below is designed to
% enforce two symmetries whose failure might disturb users:
% VALUES exactly hermitian ==> COEFFS exactly real
% VALUES exactly skew-hermitian ==> % COEFFS exactly imaginary
% This is necessary because the MATLAB FFT code does not
% exactly preserve symmetries.

% Get the length of the input:
n = size(values, 1);

% Trivial case (constant or empty):
if ( n <= 1 )
    coeffs = values; 
    return
end

% test for symmetry
vals = double([values;values(1,:)]);
isHerm = max(abs(vals-conj(vals(end:-1:1, :))),[],1) == 0;
isSkew = max(abs(vals+conj(vals(end:-1:1, :))),[],1) == 0;

% compute coefficients
coeffs = (1/n)*fftshift(fft(values, [], 1), 1);

% correct if symmetric
coeffs(:,isHerm) = real(coeffs(:,isHerm));
coeffs(:,isSkew) = 1i*imag(coeffs(:,isSkew));

% These coefficients are for interpolation defined on [0,2*pi), but we want
% to work on [-pi,pi). To fix the coefficients for this we just need to
% assign c_k = (-1)^k c_k, for k=-(N-1)/2:(N-1)/2 for N odd and 
% k = -N/2:N/2-1 for N even.
if ( mod(n, 2) ) 
    even_odd_fix = (-1).^(-(n-1)/2:(n-1)/2).';
else
    even_odd_fix = (-1).^((-n/2):(n/2-1)).';
end

coeffs = bsxfun(@times, coeffs, even_odd_fix);

end
