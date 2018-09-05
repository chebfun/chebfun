function ff = sum2(f,dim_1,dim_2)
%SUM2  Integration of a BALLFUN function over lambda and theta
%   SUM2(F,dim_1,dim_2) is the integration of the BALLFUN function f over
%   dim_1 and dim_2
%
% See also SUM, SUM3.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Sort the dimensions
dim = sort([dim_1,dim_2]);
dim_1 = dim(1); dim_2 = dim(2);

[m,n,p] = size(f);
m = m+2; p = p+2;
F = coeffs3(f,m,n,p);

if dim_1 == 2 && dim_2 == 3
    % Integrate over lambda and theta
    % Multiply by F by the measure r^2*sin(theta)

    % Extract the 0-th Fourier mode
    F = reshape(F(:,floor(n/2)+1,:), m, p);

    % Multiply f par r^2sin(theta) (= Jacobian)
    % Just do the multiplication for the 0-th Fourier mode
    Msin = trigspec.multmat(p, [.5i;0;.5i] );
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
    
elseif dim_1 == 1 &&  dim_2 == 2
    % Integrate over r and lambda
    
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
    
    % Integrate over lambda
    F = 2*pi*reshape(F(:,floor(n/2)+1,:), m, p);
    
    % Integrate over r
    A = IntChebyshev*F;
    
    % Return the chebfun function of theta on the domain
    ff = chebfun(A.', [-pi,pi],'coeffs', 'trig');

elseif dim_1 == 1 && dim_2 == 3
    % Integrate over r and theta
    
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
    
    % Permute F
    F = permute(F,[1,3,2]);
    
    % Integrate over r and theta
    A = zeros(n,1);
    
    for k = 1:n
        A(k) = IntChebyshev*F(:,:,k)*IntTheta;    
    end
        
    % Return the chebfun function of lambda on the domain
    ff = chebfun(A, [-pi,pi],'coeffs', 'trig');
else
    error('BALLFUN:sum2: Input arguments not supported'); 
end
end