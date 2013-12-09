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
%   See also COEFFS2VALS, CHEBPTS.

%   [Mathematical reference] [TODO!]

%   Copyright 2013 by The University of Oxford and The Chebfun Developers. 
%   See http://www.chebfun.org for Chebfun information.

% Get the length of the input:
n = size(values, 1);

% Trivial case (constant):
if ( n == 1 )
    coeffs = values; 
    return
end

% Flip the input data up to down. The first row of the resulting data should
% correspond to the rightmost Chebyshev point of the 1st kind in [-1 1].
values = values(end:-1:1,:);

if ( isreal(values) )
    coeffs = vals2coeffsReal(values);
elseif ( isreal(1i*values) )
    coeffs = 1i*vals2coeffsReal(imag(values));
else
    coeffs = vals2coeffsReal(real(values)) + 1i*vals2coeffsReal(imag(values));
end

end

function c = vals2coeffsReal(v)
%VALS2COEFFSREAL   Convert values at Chebyshev points to Chebyshev coefficients
%when values are real.

n = size(v,1);
m = size(v,2);
w = repmat((2/n)*exp(-1i*(0:n-1)*pi/(2*n)).',1,m);

% [TODO KX]: Confirm whether the difference between the two cases suggested
% by the remark below about "no need" is genuine.
%
% Form the vector whose data are periodic. Note that in contrast to the
% situation with 2nd kind Chebyshev points, here there is no need to mirror the
% data in order to achieve the effect of a DCT with an FFT. Even and odd cases
% are treated differently.

if ( rem(n, 2) == 0 ) % Even n case
    vv = [ v(1:2:n-1,:); v(n:-2:2,:) ];
else                  % Odd n case
    vv = [ v(1:2:n,:); v(n-1:-2:2,:) ];
end
c = real(w.*fft(vv));

% Flip back so that the trailing coefficients show up at the top rows:
c = c(end:-1:1,:);

% Halve C(N+1), i.e., the constant term in the Chebyshev series:
c(end,:) = c(end,:)/2;

end
