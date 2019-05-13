function ff = sum2(f, dims)
% SUM2 Definite integration of a BALLFUN in two variables.
%   SUM2(F, DIMS) integrates F over two of the variables r, lambda or theta 
%   where DIMS is a row vector containing two of the three indices
%   1, 2 or 3. The output is a 1D CHEBFUN in the remaining variable.
%
%   G = SUM2(F) is the same as SUM2(F, [2, 3]).
%
% See also SUM, SUM3. 

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check
if ( isempty(f) )
    ff = chebfun();
    return
end


% Default to lambda and theta directions: 
if ( nargin == 1 )
    dims = [2 3];
end

% Sort the dimensions
dims = sort(dims);
dim_1 = dims(1); dim_2 = dims(2);

% Increase the discretization by 2 in the r and theta direction
[m,n,p] = size(f);
m = m+2; p = p+2;
F = coeffs3(f,m,n,p);

if dim_1 == 2 && dim_2 == 3
    % Integrate over lambda and theta
    % Multiply by F by the measure sin(theta)

    % Extract the 0-th Fourier mode
    F = reshape(F(:,floor(n/2)+1,:), m, p);

    % Multiply f par sin(theta) (= Jacobian)
    % Just do the multiplication for the 0-th Fourier mode
    Msin = trigspec.multmat(p, [.5i;0;-.5i] );
    F = F*(Msin.');

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
    
    % Multiplication by r^2 in Chebyshev basis
    Mr2 = ultraS.multmat(m, [.5;0;.5], 0 );
    F = reshape(Mr2*reshape(F, m, []),m, n, p);
    
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
    
    % Multiplication by r^2 in Chebyshev basis
    Mr2 = ultraS.multmat(m, [.5;0;.5], 0 );
    F = reshape(Mr2*reshape(F, m, []),m, n, p);
    
    % Multiplication by theta
    Msin = trigspec.multmat(p, [.5i;0;-.5i]);
    F = reshape(reshape(F, [], p)*Msin.',m, n, p);
    
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