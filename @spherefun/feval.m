function y = feval(f, x, y, z)
%FEVAL  Evaluate a SPHEREFUN at one or more points.
%   Y = FEVAL(F, LAMBDA, THETA)  evaluates a spherefun F at (LAMBDA, THETA)
%   where LAMBDA and THETA are doubles representing the longitudinal (or 
%   azimuthal) and latitudinal (or elevation) angles.
%
%   Y = FEVAL(F, X, Y, Z)  evaluates a spherefun F at a point (X,Y,Z) in
%   Cartesian cooridnates on the surface of a sphere.
%
% See also SUBSREF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


if ( nargin == 3 )      % Spherical coordinates used.
    lambda = x;
    theta = y;
    y = feval@separableApprox(f, lambda, theta);

    if ( (size(lambda, 1) == 1) && (size(theta, 1) == 1) )
        y = y.';
    end
    
elseif ( nargin == 4 ) % Cartesian coordinates used.
    if ( isnumeric(x) && isnumeric(y) && isnumeric(z) )
        % Convert to spherical coordinates
        [lambda, theta, rad] = cart2sph(x,y,z);
        
        if ( any(rad > (1 + 1e-8)) )
            error('CHEBFUN:SPHEREFUN:FEVAL:pointsNotOnSphere',...
                ['The specified points to evaluate the function do not '...
                'lie sufficiently close to the surface of the '...
                'unit sphere.']);
        end
        % Check latitudinal coordinate system to set the elevation angle
        % appropriately.
        if iscolat(f)
            theta = pi/2 - theta;
        end
        y = feval@separableApprox(f, lambda, theta);
%         y = FastSphereEval( f, lambda, theta );
        if ( (size(lambda, 1) == 1) && (size(theta,1) == 1) )
            y = y.';
        end
        
    elseif ( strcmp(x, ':') && strcmp(y, ':') && strcmp(z, ':') )
        y = f;
        
    elseif ( strcmp(x, ':') && isnumeric(y) && isnumeric(z) )
        y = [ feval(f, sqrt(1-y.^2-z.^2), y, z) ; 
              feval(f, -sqrt(1-y.^2-z.^2), y, z) ];
        
    elseif ( strcmp(y, ':') && isnumeric(x) && isnumeric(z) )
        y = [ feval(f, x, sqrt(1-x.^2-z.^2), z); 
              feval(f, x, -sqrt(1-x.^2-z.^2), z) ];
        
    elseif ( strcmp(z, ':') && isnumeric(x) && isnumeric(y) )
        y = [ feval(f, x, y, sqrt(1-x.^2-y.^2)); 
              feval(f, x, y, -sqrt(1-x.^2-y.^2)) ];
        
    elseif ( strcmp(x, ':') && strcmp(y, ':') && isnumeric(z) )
        y = chebfun(@(t) feval(f, sqrt(1-z.^2).*cos(t), ...
                              sqrt(1-z.^2).*sin(t), z), [-pi, pi], 'trig');
        
    elseif ( strcmp(x, ':') && strcmp(z, ':') && isnumeric(y) )
        y = chebfun(@(t) feval(f, sqrt(1-y.^2).*cos(t), y, ...
                                 sqrt(1-y.^2).*sin(t)), [-pi, pi], 'trig');
        
    elseif ( strcmp(y, ':') && strcmp(z, ':') && isnumeric(x) )
        y = chebfun(@(t) feval(f, x, sqrt(1-x.^2).*cos(t), ...
                                 sqrt(1-x.^2).*sin(t)), [-pi, pi], 'trig');
        
    else
        error('CHEBFUN:SPHEREFUN:feval:argin',['Unkown input '...
            'feval(%s,%s,%s)',x,y,z]);
    end

end

end


% Fast spherefun evaluation when [M, N] = length(F) are large or there are 
% many evaluation points: 
function vals = FastSphereEval( f, lambda, theta )
% FASTSPHEREEVAL       Evaluate a spherefun using the nonuniform 2D FFT
% 
% VALS = FASTSPHEREEVAL(C, LAMBDA, THETA): 
%  
%     Inputs: 
%          F = spherefun object  
%          (LAMBDA, THETA) = Evaluation points in spherical coordinates
%     Outputs: 
%          VALS = f(LAMBDA, THETA) 
% 
% The transform costs O(m*n*(log(m) + log(n)) + N ) operations. Only use 
% use this code if [M, N] = length(F) are large or N is large. 

%%% DEVELOPERS NOTE: 
% The direct algorithm to perform this task is: 
% 
% [C, D, R] = coeffs2( f ); 
% [M, N] = size( lambda ); 
% [m, n] = length( f );
% Y = zeros(M, N); 
% mm = -floor(m/2):floor(m/2);
% nn = -floor(n/2):floor(n/2); 
% for j = 1:M
%     for k = 1:N 
%         Y(j,k) = (exp(1i*theta(j,k)*mm)*C)*D*(R.'*exp(1i*lambda(j,k)*nn'));
%     end
% end 
%
% The algorithm in this MATLAB script is based on the paper:
%
% [1] D. Ruiz--Antoln and A. Townsend, "A nonuniform fast Fourier transform
% based on low rank approximation", in preparation, 2016.
%
% Author: Alex Townsend, January 2017. 

% Convert to NUFFT2D convention (from Spherefun's convention): 
lambda = -lambda/2/pi;
theta = -theta/2/pi; 

% Working accuracy: 
tol = 10*eps; 

% Low rank approximation for the coefficients of F: 
[C, D, R] = coeffs2( f ); 
    
% Get size of evaluation points and Fourier coefficients: 
[M, N] = size( lambda ); 
[m, n] = length( f );

% Make columns: 
lambda = lambda(:); 
theta = theta(:);

% Assign each evaluation point to its closest point on a 2D equispaced
% grid:
xj = round(n*lambda)/n;
xt = mod(round(n*lambda),n)+1;
yj = round(m*theta)/m; 
yt = mod(round(m*theta),m)+1;

% Fourier modes: 
mm = -floor(m/2):floor(m/2);
nn = -floor(n/2):floor(n/2);

%%%%%%%%%%% FAST TRANSFORM %%%%%%%%%%%%%%%%
% Ay = exp(-2*pi*1i*m*(y(:)-yj(:))*mm/m);
% Find low rank approximation to Ay = U1*V1.': 
er = m*(theta-yj);
gam = norm(er, inf); 
K1 = ceil(5*gam*exp(lambertw(log(10/tol)/gam/7)));
U1 = (ChebP(K1-1,er/gam)*Bessel_cfs(K1, gam));
V1 = ChebP(K1-1, 2*mm'/m);

% Ax = exp(-2*pi*1i*n*(x(:)-xj(:))*nn/n);
% Find low rank approximation to Ay = U1*V1.': 
er = n*(lambda-xj);
gam = norm(er, inf); 
K2 = ceil(5*gam*exp(lambertw(log(10/tol)/gam/7)));
U2 = (ChebP(K2-1,er/gam)*Bessel_cfs(K2, gam));
V2 = ChebP(K2-1, 2*nn'/n);

% Business end of the transform.   (Everything above could be considered
% precomputation.) 
vals = zeros(M*N, 1);
for r = 1:K2 
    % Transform in the "row" variable: 
    Ar = ( fft(ifftshift(bsxfun(@times,(C*D*R.').',V2(:,r)),1),[],1) ).';
    for s = 1:K1
        % Transform in the "col" variable: 
        fj = fft( ifftshift(bsxfun(@times, Ar, V1(:,s)),1) );
        % Spread the love out from equispaced points to actual evaluation
        % points. Do this K1*K2 times for an accurate transform: 
        vals = vals + U1(:,s).*fj(yt(:)+m*(xt(:)-1)).*U2(:,r);
    end
end

% Reshape "vals" to the same shape as LAMBDA and THETA: 
vals = reshape( vals, M, N); 
end

function cfs = Bessel_cfs(K, gam)
% The bivarate Chebyshev coefficients for the function f(x,y) = exp(-i*x.*y)
% on the domain [-gam, gam]x[0,2*pi] are given by Lemma A.2 of Townsend's
% DPhil thesis.
arg = -gam*pi/2;
[pp,qq] = meshgrid(0:K-1);
cfs = 4*(1i).^qq.*besselj((pp+qq)/2,arg).*besselj((qq-pp)/2, arg);
cfs(2:2:end,1:2:end) = 0;
cfs(1:2:end,2:2:end) = 0;
cfs(1,:) = cfs(1,:)/2;
cfs(:,1) = cfs(:,1)/2;
end

function T = ChebP( n, x )
% Evaluate Chebyshev polynomials of degree 0,...,n at points in x. Use the
% three-term recurrence relation:
N = size(x, 1);
T = zeros(N, n+1);
T(:,1) = 1;
T(:,2) = x;
twoX = 2*x;
for k = 2:n
    T(:,k+1) = twoX.*T(:,k) - T(:,k-1);
end
end