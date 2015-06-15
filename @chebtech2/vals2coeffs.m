function coeffs = vals2coeffs(values)
%VALS2COEFFS   Convert values at Chebyshev points to Chebyshev coefficients.
%   C = VALS2COEFFS(V) returns the (N+1)x1 vector of coefficients such that F(x)
%   = C(1)*T_0(x) + C(2)*T_1(x) + C(N+1)*T_N(x) (where T_k(x) denotes the
%   k-th 1st-kind Chebyshev polynomial) interpolates the data [V(1) ; ... ;
%   V(N+1)] at Chebyshev points of the 2nd kind.
%
%   If the input V is an (N+1)xM matrix, then C = VALS2COEFFS(V) returns the
%   (N+1)xM matrix of coefficients C such that F_j(x) = C(1,j)*T_0(x) + ... +
%   C(N+1,j)*T_N(x) interpolates [V(1,j) ; ... ; V(N+1,j)] for j = 1:M.
%
% See also COEFFS2VALS, CHEBPTS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Developer Note]: This is equivalent to the Inverse Discrete Cosine Transform
% of Type I.
%
% [Mathematical reference]: Section 4.7 Mason & Handscomb, "Chebyshev
% Polynomials". Chapman & Hall/CRC (2003).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the length of the input:
n = size(values, 1);

% Trivial case (constant):
if ( n <= 1 )
    coeffs = values; 
    return
end

% Mirror the values (to fake a DCT using an FFT):
tmp = [values(n:-1:2,:) ; values(1:n-1,:)];

if ( isreal(values) )
    % Real-valued case:
    coeffs = ifft(tmp);
    coeffs = real(coeffs);
elseif ( isreal(1i*values) )
    % Imaginary-valued case:
    coeffs = ifft(imag(tmp));
    coeffs = 1i*real(coeffs);
else
    % General case:
    coeffs = ifft(tmp);
end

% Truncate:
coeffs = coeffs(1:n,:);

% Scale the interior coefficients:
coeffs(2:n-1,:) = 2*coeffs(2:n-1,:);

end
