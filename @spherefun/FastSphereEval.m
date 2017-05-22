function vals = FastSphereEval( f, lambda, theta )
% FASTSPHEREEVAL        Fast Spherefun evaluation
%
% VALS = FASTSPHEREEVAL(C, LAMBDA, THETA):
%
%     C = 2D Fourier coefficients,
%     (LAMBDA,THETA) = Evaluation points in spherical coordinates
%     VALS = Values of Fourier series with coefficients C at (LAMBDA,THETA)
%
%  The transform costs O(K2*m*n*log(m*n) + K1*K2*N ) operations.
%
%  This is a fast transform analogue to the direct algorithm:
% [M, N] = size( x );
% [m, n] = size( c );
% Y = zeros(M, N);
% mm = -floor(m/2):floor(m/2);
% nn = -floor(n/2):floor(n/2);
% for j = 1:M
%     for k = 1:N
%         Y(j,k) = exp(1i*theta(j,k)*mm)*c*exp(1i*lambda(j,k)*nn');
%     end
% end
%
% The algorithm in this MATLAB script is based on the paper:
%
% [1] D. Ruiz--Antoln and A. Townsend, "A nonuniform fast Fourier transform
% based on low rank approximation", in preparation, 2016.
%
% This paper related the NUFFT to a FFT by low rank approximation.
%
% Author: Alex Townsend, January 2017.

% Convert to NUFFT2D convention (from Spherefun's convention):
lambda = -lambda/2/pi;
theta = -theta/2/pi;

% Working accuracy:
tol = 10*eps;

% Get size of evaluation points and Fourier coefficients:
[M, N] = size( lambda );
[n, m] = length( f );
m = m + (1-mod(m,2));
n = n + (1-mod(n,2));

% Low rank approximation for the coefficients of F:
[C, D, R] = coeffs2( f, n, m);
rk = size(D,1);

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
% K1 = ceil(5*gam*exp(lambertw(log(10/tol)/gam/7)));
K1 = 15;
U1 = (ChebP(K1-1,er/gam)*Bessel_cfs(K1, gam));
V1 = ChebP(K1-1, 2*mm'/m);

% Ax = exp(-2*pi*1i*n*(x(:)-xj(:))*nn/n);
% Find low rank approximation to Ay = U1*V1.':
er = n*(lambda-xj);
gam = norm(er, inf);
% K2 = ceil(5*gam*exp(lambertw(log(10/tol)/gam/7)));
K2 = 15;
U2 = (ChebP(K2-1,er/gam)*Bessel_cfs(K2, gam)).';
V2 = ChebP(K2-1, 2*nn'/n);

% Business end of the transform.   (Everything above could be considered
% precomputation.)
FFT_cols = zeros(m,rk,K1);
for s = 1:K1
    % Transform in the "col" variable:
    FFT_cols(:,:,s) = fft( ifftshift(bsxfun(@times, C, V1(:,s)),1) );
end
FFT_rows = zeros(rk,n,K2);
for r = 1:K2
    % Transform in the "row" variable:
    FFT_rows(:,:,r) = D*fft(ifftshift(bsxfun(@times,R,V2(:,r)).',2),[],2);
end

XX = permute(FFT_rows,[1 3 2]);
YY = permute(FFT_cols,[3 2 1]);

% Spread the love out from equispaced points to actual evaluation
% points. Do this K1*K2 times for an accurate transform:
[ii, jj] = ind2sub( [m,n], yt(:)+m*(xt(:)-1) );
% Aij = zeros(K1, K2, numel(ii) );

% Convert to cell for speed:
XXX = cell(size(XX,3),1);
for k = 1:length(XXX)
    XXX{k} = XX(:,:,k);
end
YYY = cell(size(YY,3),1);
for k = 1:length(YYY)
    YYY{k} = YY(:,:,k);
end
% Do only unique multiplications:
[c, ~, ic2] = unique([ii jj], 'rows');
temp = cell(size(c,1),1);
for idx = 1:size(c,1)
    temp{idx} = (YYY{c(idx,1)}*XXX{c(idx,2)});
end
% Recover A:
vals = zeros(numel(ii),1);
for idx = 1:numel(ii)
    vals(idx) = U1(idx,:)*(temp{ic2(idx)}*U2(:,idx));
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