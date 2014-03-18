function values = coeffs2vals(coeffs)
%COEFFS2VALS   Convert Fourier coefficients to values at equally spaced
%points between [-pi pi).
%   V = COEFFS2VALS(C) returns the values of the trignometric polynomial 
%   as follows:
%   If N is odd
%   V(i,1) = P(x_i) =
%   C(1,1)*z^{-(N-1)/2} + C(2,1)*z^{-(N-1)/2+1} + ... + C(N,1)*z^{(N-1)/2}, where the x_i are
%
%   [FIXME] 2nd-kind Chebyshev nodes.
%
%   If the input C is an (N+1)xM matrix then V = COEFFS2VALS(C) returns the
%   (N+1)xM matrix of values V such that V(i,j) = P_j(x_i) = C(1,j)*T_{N-1}(x_i)
%   + C(2,j)*T_{N-2}(x_i) + ... + C(N,j)
%
% See also VALS2COEFFS, FOURIERPTS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Get the length of the input:
n = size(coeffs, 1);

% Trivial case (constant):
if ( n == 1 )
    values = coeffs; 
    return
end

% Shift the coefficients properly
values = ifft(ifftshift(flipud(n*coeffs),1),[],1);

end
