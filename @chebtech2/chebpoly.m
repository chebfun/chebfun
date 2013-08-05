function coeffs = chebpoly(values)
%CHEBPOLY   Convert values at Chebyshev points to Chebyshev coefficients.
%   C = CHEBPOLY(V) returns the (N+1)x1 vector of coefficients such that F(x) =
%   C(1)*T_N(x) + ... + C(N)*T_1(x) + C(N+1)*T_0(x) (where T_k(x) denotes the
%   k-th 1st-kind Chebyshev polynomial) interpolates the data [V(1) ; ... ;
%   V(N+1)] at Chebyshev points of the 2nd kind.
%
%   If the input V is an (N+1)xM matrix, then C = CHEBPOLY(V) returns the
%   (N+1)xM matrix of coefficients C such that F_j(x) = C(1,j)*T_N(x) + ... +
%   C(N,j)*T_1(x) + C(N+1,j)*T_0(x) interpolates [V(1,j) ; ... ; V(N+1,j)] for j
%   = 1:M.
%
%   See also CHEBPOLYVAL, CHEBPTS.

%   Copyright 2013 by The University of Oxford and The Chebfun Developers. 
%   See http://www.chebfun.org for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   [Mathematical reference]: Section 4.7 Mason & Handscomb, "Chebyshev
%   Polynomials". Chapman & Hall/CRC (2003).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the length of the input:
n = size(values, 1);

% Trivial case (constant):
if ( n == 1 )
    coeffs = values; 
    return
end

% Mirror the values (to fake a DCT using an FFT):
tmp = [values(end:-1:2,:) ; values(1:end-1,:)];

if ( isreal(values) )
    % Real-valued case:
    coeffs = fft(tmp)/(n-1);
    coeffs = real(coeffs);
elseif ( isreal(1i*values) )
    % Imaginary-valued case:
    coeffs = fft(imag(tmp))/(n-1);
    coeffs = 1i*real(coeffs);
else
    % General case:
    coeffs = fft(tmp)/(n-1);
end

% Truncate and flip the order:
coeffs = coeffs(n:-1:1,:);

% Scale the first and final coefficient:
coeffs([1,end],:) = .5*coeffs([1,end],:);

end
