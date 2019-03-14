function g = sum(f, dim)
% SUM Definite integration of a BALLFUN.
%   SUM(F, DIM) where DIM is 1, 2 or 3 integrates only over r (radial direction), 
%   lambda (azimuthal direction) or theta (polar direction) respectively and 
%   and returns as its output a spherefun if DIM is 1 or a diskfun otherwise.
%
% See also SUM2, SUM3. 

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Default direction: 
if ( nargin == 1 )
    dim = 1;
end

if isempty( f )
    if dim == 1
        g = spherefun();
    elseif (dim == 2) || (dim == 3)
        g = diskfun();
    % Return an error
    else
        error('BALLFUN:sum: Input arguments not supported'); 
    end
    return
end

% Increase the discretization by 2 in the r and theta direction
[m,n,p] = size(f);
m = m+2; p = p+2;
F = coeffs3(f,m,n,p);

% Compute the integral over r of r^2*f
if dim == 1
    % Multiplication by r^2 in Chebyshev basis
    Mr2 = ultraS.multmat(m, [.5;0;.5], 0);
    F = reshape(Mr2*reshape(F, m, []),m, n, p);
    
    % Coefficients of integration of the Chebyshev polynomials
    K = zeros(1,m);
    for i = 0:m-1
        if mod(i,4)==0
            K(i+1) = -1/(i^2-1);
        elseif mod(i,4)==1
            K(i+1) = 1/(i+1);
        elseif mod(i,4)==2
            K(i+1) = -1/(i^2-1);
        else
            K(i+1) = -1/(i-1);
        end
    end

    % Matrix of coefficients of Fourier-Fourier after integration over r
    A = zeros(n,p);
    for i = 1:m
       A = A+K(i).*reshape(F(i,:,:),n, p);
    end

    % Return the spherefun function of lambda,theta; coeffs2spherefun takes a
    % matrix of coefficients of theta and lambda instead of lambda, theta
    g = spherefun.coeffs2spherefun(A.');
 
% Compute the integral over lambda of f
elseif dim ==2
    % Integration for the Fourier modes in lambda
    Listp = (1:n) - floor(n/2) - 1;
    C = -1i*((-1).^Listp-1)./Listp;
    C(floor(n/2)+1) = 2*pi;
    
    % Matrix of coefficients of Chebyshev-Fourier after integration over
    % lambda
    A = zeros(m,p);
    for i = 1:n
       A = A + C(i)*reshape(F(:,i,:),m, p); 
    end
    
    % Return the diskfun function of r, theta
    g = diskfun.coeffs2diskfun(real(A));
    
% Compute the integral over theta of f*sin(theta)
elseif dim == 3
    
    % Multiplication by theta
    Msin = trigspec.multmat(p, [.5i;0;-.5i]);
    F = reshape(reshape(F, [], p)*Msin.',m, n, p);
    
    % Integration for the Fourier modes in theta
    Listp = (1:p) - floor(p/2) - 1;
    C = -1i*((-1).^Listp-1)./Listp;
    C(floor(p/2)+1) = pi;
    % Matrix of coefficients of Chebyshev-Fourier after integration over
    % lambda
    A = zeros(m,n);
    for i = 1:p
       A = A + C(i)*F(:,:,i);
    end
    
    % Return the diskfun function of r, lambda
    g = diskfun.coeffs2diskfun(real(A));
        
% Return an error
else
    error('BALLFUN:sum: Input arguments not supported'); 
end
end
