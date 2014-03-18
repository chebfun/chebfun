function coeffs = vals2coeffs(values)
%VALS2COEFFS   Convert values at equally spaced points between [-pi pi).
%   TODO: NEEDS UPDATING
%   C = VALS2COEFFS(V) returns the (N+1)x1 vector of coefficients such that F(x) =
%   C(1)*T_N(x) + ... + C(N)*T_1(x) + C(N+1)*T_0(x) (where T_k(x) denotes the
%   k-th 1st-kind Chebyshev polynomial) interpolates the data [V(1) ; ... ;
%   V(N+1)] at Chebyshev points of the 2nd kind.
%
%   If the input V is an (N+1)xM matrix, then C = VALS2COEFFS(V) returns the
%   (N+1)xM matrix of coefficients C such that F_j(x) = C(1,j)*T_N(x) + ... +
%   C(N,j)*T_1(x) + C(N+1,j)*T_0(x) interpolates [V(1,j) ; ... ; V(N+1,j)] for j
%   = 1:M.
%
% See also COEFFS2VALS, CHEBPTS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Get the length of the input:
n = size(values, 1);

% Trivial case (constant):
if ( n <= 1 )
    coeffs = values; 
    return
end

coeffs = (1/n)*flipud(fftshift(fft(values,[],1),1));

end
