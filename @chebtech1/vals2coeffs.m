function coeffs = vals2coeffs(values)
%VALS2COEFFS   Convert values at Chebyshev points to Chebyshev coefficients.
%   C = VALS2COEFFS(V) returns the (N+1)x1 vector of coefficients such that
%   F(x) = C(1)*T_N(x) + ... + C(N)*T_1(x) + C(N+1)*T_0(x) (where T_k(x)
%   denotes the k-th 1st-kind Chebyshev polynomial) interpolates the data
%   [V(1) ; ... ; V(N+1)] at Chebyshev points of the 1st kind. 
%
%   If the input V is an (N+1)xM matrix, then C = VALS2COEFFS(V) returns the
%   (N+1)xM matrix of coefficients C such that F_j(x) = C(1,j)*T_N(x) + ... 
%   + C(N,j)*T_1(x) + C(N+1)*T_0(x) interpolates [V(1,j) ; ... ; V(N+1,j)]
%   for j = 1:M.
%
% See also COEFFS2VALS, CHEBPTS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Developer Note]: This is equivalent to Discrete Cosine Transform of Type II.
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

% Computing the weight vector often accounts for at least half the cost of this
% transformation. Given that (a) the weight vector depends only on the length of
% the coefficients and not the coefficients themselves and (b) that we often
% perform repeated transforms of the same length, we store w persistently.
persistent w
if ( size(w, 1) ~= n )
    % Pre-compute the weight vector:
    w = 2*exp(1i*(0:n-1)*pi/(2*n)).';
end

% Mirror the values for FFT:
tmp = [values(n:-1:1, :) ; values];
coeffs = ifft(tmp);

% Truncate, flip the order, and multiply the weight vector:
coeffs = bsxfun(@times, w, coeffs(1:n, :));

% Scale the coefficient for the constant term:
coeffs(1,:) = coeffs(1,:)/2;

% Post-process:
if ( isreal(values) )  
    % Real-valued case:
    coeffs = real(coeffs);
elseif ( isreal(1i*values) )  
    % Imaginary-valued case:
    coeffs = 1i*imag(coeffs);
end

end
