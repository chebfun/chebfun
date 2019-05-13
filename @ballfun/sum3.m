function I = sum3(f)
%SUM3   Triple integral of a BALLFUN over its domain.
%   I = SUM3(F) returns the double definite integral of a BALLFUN.
%
% See also SUM, SUM2.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check 
if ( isempty(f) )
    I = [];
    return
end

% Second argument: integral over the sphere of radius -1
[m,n,p] = size(f);
m = m+2; p = p+2;
F = coeffs3(f,m,n,p);

% Extract the 0-th Fourier mode
F = reshape(F(:,floor(n/2)+1,:), m, p);

% Increase the discretization by 2 in the r and theta direction
F = [zeros(m,1),F,zeros(m,1);zeros(2,p+2)];
m = m+2;
p = p+2;

% Multiply f par r^2sin(theta) (= Jacobian)
Msin = trigspec.multmat(p, [0.5i;0;-0.5i] );
Mr2 = ultraS.multmat(m, [0.5; 0; 0.5], 0 );
F = Mr2*F*(Msin.');

% Coefficients of integration between 0 and 1 of the chebyshev polynomials
IntChebyshev = zeros(1,m);
for i = 0:m-1
    if mod(i,4)==0
        IntChebyshev(i+1) = -1/(i^2-1);
    elseif mod(i,4)==1
        IntChebyshev(i+1) = 1/(i+1);
    elseif mod(i,4)==2
        IntChebyshev(i+1) = -1/(i^2-1);
    else
        IntChebyshev(i+1) = -1/(i-1);
    end
end

% Coefficients of integration between 0 and pi of the theta Fourier
% function
Listp = (1:p).' - floor(p/2)-1;
IntTheta = -1i*((-1).^Listp-1)./Listp;
IntTheta(floor(p/2)+1) = pi;

% Integrate over lambda
IntTheta = 2*pi*IntTheta;

% Integrate over the sphere of radius -1
if nargin>1
    IntChebyshev = (-1).^(0:m-1).*IntChebyshev;
end

% Return the integral of f over the ballfun
I = IntChebyshev*F*IntTheta;

% Return real value if the function is real
if f.isReal
   I = real(I); 
end
end
