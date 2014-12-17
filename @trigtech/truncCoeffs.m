function c = truncCoeffs(f, n)
%TRUNCCOEFFS   coefficients of degree N least square approximation of F.
%   C = TRUNCCOEFFS(F, N) returns the middle 2N+1 coefficients of F.

c = f.coeffs;
N = length(c);

% If required degree n is larger than the existing coefficients, return
% all:
if ( N <= 2*n + 1 )
    return
end

% Discard the onese not needed:
if ( rem(N,2) == 1)
    k = (N+1)/2;
else
    k = N/2+1;
end
c = [c(k-n:k-1); c(k:k+n)];

end

