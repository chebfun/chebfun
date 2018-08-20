function ff = sum2(f)
% SUM2 Integration of a BALLFUN function over lambda and theta
%   SUM2(f) is the integration of the BALLFUN function f over lambda 
%   and theta

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[m,n,p] = size(f);

F = f.coeffs;

% Extract the 0-th Fourier mode
F = reshape(F(:,floor(n/2)+1,:), m, p);

% Multiply f par r^2sin(theta) (= Jacobian)
% Just do the multiplication for the 0-th Fourier mode
trig1 = trigtech( @(t) sin(pi*t));
Msin = trigspec.multmat(p, trig1.coeffs );
Mr2 = ultraS.multmat(m, [.5;0;.5], 0 );
F = Mr2*F*(Msin.');

% Integration for the Fourier modes in theta
Listp = (1:p) - floor(p/2) - 1;
C = -1i*((-1).^Listp-1)./Listp;
C(floor(p/2)+1) = pi;

C=2*pi*C.';
% Vector of coefficients of chebyshev after integration over lambda and
% theta
A = F*C;

% Return the chebfun function of r after the multiplication by r^2
% (measure)
ff = chebfun(A, 'coeffs');
end
