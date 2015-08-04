function f = constructorTurbo(f, op, pref)
%   F = CONSTRUCTORTURBO(F, OP, PREF) takes the CHEBTECH F that was constructed
%   from the function handle OP using the preferences in PREF and tries to use
%   contour integrals to compute some additional higher-order coefficients to
%   high accuracy.  These coefficients are returned in the appropriate
%   positions in F.COEFFS.
%
%   This function is only meant to be called by the CHEBTECH constructor.
%
% References:
%
%   [1] Bornemann, F.  Accuracy and stability of computing high-order
%         derivatives of analytic functions by Cauchy integrals.  Found.
%         Comput. Math. 11 (2011), pp. 1-63.
%
%   [2] Wang, H. and Huybrechs, D.  Fast and accurate computation of Jacobi
%         expansion coefficients of analytic functions.  Technical Report
%         TW-645, K.U. Leuven, 2015.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% How many coefficients do we want to compute?
if ( ~isnan(pref.fixedLength) )
    n = pref.fixedLength;
else
    n = 2*length(f);
end

% Pick an ellipse over which to integrate.
rhoCheb = exp(abs(log(eps())) / length(f));
rho = rhoCheb^(2/3);

% Do the contour integrals.
c = chebCoeffsTurbo(op, rho, n);

% Assign the new coefficients.
if ( isreal(f) )
    f.coeffs = real(c);
elseif ( isreal(1i*f) )
    f.coeffs = imag(c);
else
    f.coeffs = c;
end

end

function c = chebCoeffsTurbo(f, rho, n)
%CHEBCOEFFSTURBO   Compute Chebyshev coefficients using contour integrals.
%   C = CHEBCOEFFSTURBO(F, RHO, N) computes the first N Chebyshev coefficients
%   of the holomorphic function represented by the function handle F using
%   Cauchy integrals over the Bernstein ellipse of parameter RHO.  F must be
%   vectorized and able to accept complex inputs.

K = 4*n;                                    % Number of quadrature nodes.
g = @(z) f((rho*z + 1./(rho*z))/2);         % Remap ellipse to unit circle.
z = exp(2*pi*1i*(0:1:(K - 1)).'/K);         % Roots of unity.

% Compute integrals with trap. rule and rescale to get coefficients.
c = bsxfun(@rdivide, fft(g(z))/K, rho.^(0:1:(K - 1)).');
c = [c(1,:) ; 2*c(2:n,:)];

end
