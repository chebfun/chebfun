function f = nufft( c, x, omega, tol)
% CHEBFUN.NUFFT   Nonuniform fast Fourier transform
%
% F = CHEBFUN.NUFFT( C ) is the same as fft( C ). C must be a column
% vector.
% 
% F = CHEBFUN.NUFFT( C, X ) is a nonuniform fast Fourier transform of type
% 1,
% which computes the following sums in quasi-optimal complexity:
% 
%       F_j = \sum_{k=0}^{N-1} C(K)*exp(-2*pi*1i*X(j)*k/N), 0<=j<=N-1.
%
% C and X must be column vectors of the same length. 
%
% F = CHEBFUN.NUFFT( C, X, 1 ) is the same as CHEBFUN.NUFFT( C, X )
%
% F = CHEBFUN.NUFFT( C, OMEGA, 2 ) is a nonuniform fast Fourier transform 
% of type 2, which computes the following sums in quasi-optimal complexity:
% 
%       F_j = \sum_{k=0}^{N-1} C(K)*exp(-2*pi*1i*j*OMEGA(k)/N), 0<=j<=N-1.
% 
% F = CHEBFUN.NUFFT( C, X, OMEGA ) is a nonuniform fast Fourier transform
% of type 3, which computes the following sums in quasi-optimal complexity:
% 
%       F_j = \sum_{k=0}^{N-1} C(K)*exp(-2*pi*1i*X(j)*OMEGA(k)/N), 0<=j<=N-1.
%
% F = CHEBFUN.NUFFT( C, X, 1, TOL), F = CHEBFUN.NUFFT(C, OMEGA, 2, TOL),
% and F = CHEBFUN.NUFFT( C, X, OMEGA, TOL ) are the same as above but with
% a tolerance of TOL. By default, TOL = eps. 
% 
% See also chebfun.ndct, chebfun.dct, chebfun.dst.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%% DEVELOPER'S NOTE %%%
% The algorithm in this MATLAB script is based on the paper: 
% [1] D. Ruiz--Antoln and A. Townsend, "A nonuniform fast Fourier transform
% based on low rank approximation", in preparation, 2016. 

if ( nargin == 1 ) 
    f = fft( c ); 
elseif ( nargin == 2 ) 
    % default to type 1 nufft
    f = nufft1( c, x, tol );
elseif ( nargin == 3 ) 
    tol = eps;
    type = omega;
    if numel(type) 
        if type == 1
            f = nufft1( c, x, tol);
        elseif type == 2
            f = nufft2( c, x, tol);
        else 
            error('CHEBFUN::NUFFT::TYPE','Unrecognised NUFFT type.');
        end
    else
        error('CHEBFUN::NUFFT:SYNTAX','Unrecognised number of arguments to NUFFT.')
    end
elseif ( nargin == 4 ) 
    
end
end
        
function f = nufft1( c, x, tol )
% NUFFT-1:
[~, t, gam] = FindAlgorithmicParameters( x );
K = ceil(5*gam*exp(lambertw(log(10/tol)/gam/7)));
[u, v, scl] = constructAK( x, K );

g = fft( repmat(c,1,K).*v );
f = scl.*sum(u.*g(t,:),2);
end

function [s, t, gam] = FindAlgorithmicParameters( x )

N = size(x, 1);
s = round(N*x);
t = mod(s, N) + 1;
gam = norm( N*x - s, inf);
end

function [u, v, scl] = constructAK( x, K )
% Construct a low rank approximation to 
%
%     A_{jk} = exp(-2*pi*1i*(x_j-j/N)*k), 0<=j,k<=N-1,
% 
% where |x_j-j/N|<= gam <=1/2.  See [1]. 

N = size(x, 1);
[s, ~, gam] = FindAlgorithmicParameters( x );
er = N*x - s;
scl = exp(-1i*pi*er);
u = ChebP(K-1,er/gam)*Bessel_cfs(K, gam); 
v = ChebP(K-1, 2*(0:N-1)'/N - 1);
end

function cfs = Bessel_cfs(K, gam)
% The bivarate Chebyshev coefficients for the function f(x,y) = exp(-i*x.*y)
% on the domain [-gam, gam]x[0,2*pi] are given by Lemma A.2 of Townsend's
% DPhil thesis.
arg = -gam*pi/2;
[pp,qq] = meshgrid(0:K-1); 
cfs = 4*(1i).^qq.*besselj((pp+qq)/2,arg).*besselj((qq-pp)/2,arg);
cfs(2:2:end) = 0;
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