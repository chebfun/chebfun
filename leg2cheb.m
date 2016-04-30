function c_cheb = leg2cheb(c_leg, varargin)
%LEG2CHEB   Convert Legendre coefficients to Chebyshev coefficients.
%   C_CHEB = LEG2CHEB(C_LEG) converts the vector C_LEG of Legendre coefficients
%   to a vector C_CHEB of Chebyshev coefficients such that
%       C_CHEB(1)*T0 + ... + C_CHEB(N)*T{N-1} = ...
%           C_LEG(N)*P0 + ... + C_LEG(1)*P{N-1},
%   where P{k} is the degree k Legendre polynomial normalized so that max(|P{k}|
%   = 1.
%
%   C_CHEB = LEG2CHEB(C_LEG, 'norm') is as above, but with the Legendre
%   polynomials normalized to be orthonormal.
%
%   C = LEG2CHEB(C_LEG, 'trans') returns the `transpose' of the LEG2CHEB
%   operator applied to C_LEG. That is, if C_CHEB = B*C_LEG, then C = B'*C_LEG.
%
%   If C_LEG is a matrix then the LEG2CHEB operation is applied to each column.
%
%   For N > 513 the algorithm used is the one described in [1].
%
%   References:
%     [1] A. Townsend, M. Webb, and S. Olver, "Fast polynomial transforms based 
%         on Toeplitz and Hankel matrices", submitted, 2016.
%
% See also CHEB2LEG, CHEB2JAC, JAC2CHEB, JAC2JAC.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[N, n] = size(c_leg); % Number of columns.

normalize = false;    % Default - no normalize.
if ( any(strncmpi(varargin, 'normalize', 4)) )
    c_leg = bsxfun(@times, c_leg, sqrt((0:N-1)'+1/2) );
end

trans = 0; 
if ( any(strncmpi(varargin, 'trans', 4)) )
    trans = true;
end

if ( N < 2 ) 
    % Trivial case:
    c_cheb = c_leg; 
elseif ( N <= 512 ) 
    % Use direct approach: 
    L = legvandermonde(N-1, cos(pi*(0:(N-1))'/(N-1)));
    if ( ~trans )
        c_cheb = idct1(L*c_leg);
    else
        c_cheb = L.'*idct1(c_leg);
    end
else
    c_cheb = leg2cheb_fast(c_leg, trans);
end

end

function c_cheb = leg2cheb_fast(c_leg, trans)
% BRIEF IDEA: 
%  Let A be the upper-triangular conversion matrix. We observe that A can be
%  decomposed as A = D1(T.*H)D2, where D1 and D2 are diagonal, T is Toeplitz,
%  and H is a real, symmetric, positive definite Hankel matrix. The Hankel part
%  can be approximated, up to an error of tol, by a rank O( log N log(1/tol) )
%  matrix. A low rank approximation is constructed via pivoted Cholesky
%  factorization.
%
% For more information on the algorithm, see Section 5.3 of [1].

[N, n] = size(c_leg); % Number of columns.

%%%%%%%%%%%%%%%%%%%%%% Initialise  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate the symbol of the Hankel part of M:

% This for-loop is a faster and more accurate way of doing:
%   Lambda = @(z) exp(gammaln(z+1/2) - gammaln(z+1));
%   vals = Lambda( (0:2*N-1)'/2 );)
vals = [sqrt(pi) ; 2/sqrt(pi) ; zeros(2*N-2,1)];
for i = 2:2:2*(N-1)
    vals(i+1) = vals(i-1)*(1-1/i);
    vals(i+2) = vals(i)*(1-1/(i+1));
end

%%%%%%%%%%%%%%%%%%  Pivoted Cholesky algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate quasi-SVD of hankel part by using Cholesky factorization.
% Find the numerical rank and pivot locations. This is equivalent to
% Cholesky factorization on the matrix A with diagonal pivoting, except
% here only the diagonal of the matrix is updated.
d = vals(1:2:2*N);    % Diagonal of Hankel matrix.
pivotValues = [];     % Store Cholesky pivots.
C = [];               % Store Cholesky columns.
tol = 1e-14;          % Tolerance of low rank approx.
k = 0;
[mx, idx] = max( d ); % Max on diagonal, searching for first Cholesky pivot.
while ( mx > tol )
    
    newCol = vals( idx:idx+N-1 );
    if ( size(C, 2) > 0)
        newCol = newCol - C*(C(idx,:).' .* pivotValues);
    end

    pivotValues = [pivotValues ; 1./mx]; %#ok<AGROW> % Append pivtoValues.
    C = [C, newCol]; %#ok<AGROW>                     % Append newCol to C.
    d = d - newCol.^2 ./ mx;                         % Update diagonal.
    [mx, idx] = max(d);                              % Find next pivot.
    
end
sz = size(C, 2);                                     % Numerical rank of H.
C = C * spdiags(sqrt(pivotValues), 0, sz, sz);       % Share out scaling.

%%%%%%%%%%%%%%%%%%  Multiply  D1(T.*H)D2  %%%%%%%%%%%%%%%%%%%%%%%%%%
% Upper-triangular Toeplitz matrix in A = D1(T.*H)D2:
T_row1 = vals(1:N);  
T_row1(2:2:end) = 0;
Z = zeros(N, 1);

% Fast Toeplitz matrix multiply. (This is the optimized since this is the
% majority of the cost of the code.)

% Toeplitz symbol: 
if ( ~trans ) 
    % Upper-triangular toeplitz: 
    a = fft( [T_row1(1) ; Z ; T_row1(N:-1:2)] );
else
    % Lower-triangular toeplitz: 
    a = fft( [T_row1(1:N) ; Z] );
    % Diagonal-scaling is on the right if trans = 0. 
    c_leg = 2/pi*c_leg; 
    c_leg(1,:) = c_leg(1,:)/2; 
end

if ( n == 1 )     % Column input
    tmp1 = bsxfun(@times, C, c_leg);
    f1 = fft( tmp1, 2*N );
    tmp2 = bsxfun(@times, f1, a );
    b = ifft( tmp2 );
    c_cheb = sum(b(1:N,:).*C, 2);
else              % Matrix input
    c_cheb = zeros(N, n);
    for k = 1:n
        tmp1 = bsxfun(@times, C, c_leg(:,k));
        f1 = fft( tmp1, 2*N, 1 );
        tmp2 = bsxfun(@times, f1, a);
        b = ifft( tmp2, [], 1 );
        c_cheb(:,k) = sum(b(1:N,:).*C, 2);
    end
end

% Diagonal-scaling is on the left if trans = 0. 
if ( ~trans ) 
    c_cheb = 2/pi*c_cheb; 
    c_cheb(1,:) = c_cheb(1,:)/2; 
end

end

function L = legvandermonde(N, x)
% Legendre-Chebyshev Vandemonde matrix:
Pm2 = 1; Pm1 = x;                           % Initialise.
L = zeros(length(x), N+1);                  % Vandermonde matrix. 
L(:,1:2) = [1+0*x, x];                      % P_0 and P_1.     
for n = 1:N-1                               % Recurrence relation:
    P = (2-1/(n+1))*Pm1.*x - (1-1/(n+1))*Pm2;  
    Pm2 = Pm1; Pm1 = P; 
    L(:,2+n) = P;
end
end

function c = idct1(v)
%IDCT1   Convert values on a Cheb grid to Cheb coefficients (inverse DCT1).
% IDCT1(V) returns T(X)\V, where X = cos(pi*(0:N)/N), T(X) = [T_0, T_1, ...,
% T_N](X) (where T_k is the kth 1st-kind Chebyshev polynomial), and N =
% length(V) - 1.

c = chebfun.idct(v, 1);                     % IDCT-I
c([1,end],:) = .5*c([1,end],:);             % Scale.

end
