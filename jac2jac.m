function c_jac = jac2jac( c_jac, alpha, beta, gam, delta )
%JAC2JAC   Convert Jacobi (A,B) coefficients to Jacobi (G,D) coefficients.
%   C_JAC = JAC2JAC(C_JAC, A, B, G, D) converts the vector C_JAC of Jacobi
%   (A,B) coefficients to a vector of Jacobi (G,D) coefficients such that
%    C_JAC(1)*P_0^{(A,B)}(x) + ... + C_JAC(N)*P_{N-1}^{(A,B)}(x) = ...
%           C_JAC(1)*P_0^{(G,D)}(x) + ... + C_JAC(N)*P{N-1}^{(G,D)}(x),
%   where P_k^{(A,B)} is the degree k Jacobi polynomial corresponding to
%   the weight function w(x) = (1-X)^A * (1+X)^B.
%
% See also JAC2CHEB, CHEB2JAC.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%% For developers  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For more information on the algorithm, see Section 5.3 of
%
% A. Townsend, M. Webb, and S. Olver, "Fast polynomial transforms based on
% Toeplitz and Hankel matrices", submitted, 2016.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vectorize the code here:
% TODO: Improve the vectorized implementation.
NumCols = size(c_jac,2);
if ( NumCols > 1 )
    for k = 1 : NumCols
        c_jac(:,k) = jac2jac( c_jac(:,k), alpha, beta, gam, delta );
        return
    end
end

% Move (alpha,beta) to (A,B) so that |A-alpha|<1 and |B-beta|<1
while ( alpha <= gam - 1 )
    c_jac = RightJacobiConversion(c_jac, alpha, beta);
    alpha = alpha + 1;
end

while ( alpha >= gam + 1 )
    c_jac = LeftJacobiConversion(c_jac, alpha-1, beta);
    alpha = alpha - 1;
end

while ( beta <= delta - 1 )
    c_jac = UpJacobiConversion(c_jac, alpha, beta);
    beta = beta + 1;
end

while ( beta >= delta + 1 )
    c_jac = DownJacobiConversion(c_jac, alpha, beta-1);
    beta = beta - 1;
end

% Now take (alpha,beta) to (gamma,delta), where |alpha-gamma|<1 and
% |beta-delta|<1:
if ( abs( alpha - gam ) > 1e-15 )
    c_jac = JacobiFractionalConversion( c_jac, alpha, beta, gam );
end

if ( abs( beta - delta ) > 1e-15 )
    % Use reflection formula for Jacobi polynomials:
    c_jac(2:2:end) = -c_jac(2:2:end);
    c_jac = JacobiFractionalConversion( c_jac, beta, gam, delta );
    c_jac(2:2:end) = -c_jac(2:2:end);
end
end

function c_jac = JacobiFractionalConversion( v, alpha, beta, gam )
%JACOBIFRACTIONALCONVERSION     Convert Jacobi (alpha,beta) -> (gam,beta).
%
% W = JACOBIFRACTIONALCONVERSION(V, A, B, G) converts the coefficients V in
% a P^(A,B) Jacobi expansion to coefficients in a (G,B) Jacobi expansion.
% Here, we assume that |A-G|<1.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% IDEA:
% Let A be the upper-triangular conversion matrix. We observe that A can be
% decomposed as A = D1(T.*H)D2, where D1 and D2 are diagonal, T is
% Toeplitz, and H is a real, symmetric, positive definite Hankel matrix.
% The Hankel part can be approximated, up to an error of tol, by a
% rank O( log N log(1/tol) ) matrix. A low rank approximation is
% constructed via pivoted Cholesky factorization.

% Go to Chebyshev if no fourth argument:
if ( nargin < 4 )
    gam = -1/2;
end

N = size(v, 1);  % Size of coefficients.

% Functions to define Toeplitz and Hankel part of A:
Lambda1 = @(z) exp( gammaln( z+alpha+beta+1 ) - gammaln(z+gam+beta+2) );
Lambda2 = @(z) exp( gammaln( z+alpha-gam) - gammaln(z+1) );
Lambda3 = @(z) exp( gammaln( z+gam+beta+1) - gammaln(z+beta+1) );
Lambda4 = @(z) exp( gammaln( z+beta+1 ) - gammaln( z + alpha +beta + 1) );

% Diagonal matrices in A = D1(T.*H)D2:
D1 = spdiags( ((2*(0:N-1)+gam+beta+1).*Lambda3([1 1:N-1]))',0,N,N ); D1(1,1) = 1;
D2 = 1./gamma(alpha-gam)*spdiags( Lambda4([1 1:N-1])', 0, N, N );
D2(1,1) = 0; 

% Symbol of the Hankel part:
vals = Lambda1( [1 1:2*N-1] )';
vals(1) = 0; 
% Note that we would then usually do the following, but it is too slow:
% H = hankel( vals(1:N), vals(N:2*N-1) );

% PIVOTED CHOLESKY ALGORITHM ON H:
% Equivalent to GE with complete pivoting on H, but much(!) faster.
% Only the diagonal of H and the pivoting columns are updated.

d = vals(1:2:end);  % Diagonal of H, equivalent to diag( H )
pivotValues = [];   % store Cholesky pivots.
C = [];             % store Cholesky columns.
tol = 1e-14*log(N); % Tolerance of low rank approx.
k = 0;
[mx, idx] = max( d ); % max on diagonal, searching for first Cholesky pivot.
while ( mx > tol )
    
    k = k + 1;
    
    pivotValues = [pivotValues ; 1./mx]; %#ok<AGROW> 
    
    % Extract column selected by pivoting.
    newCol = vals(idx:idx+N-1); % Equivalent to H(:,idx)
    
    for j = 1 : k-1
        newCol = newCol - C(:, j) * C(idx, j) * pivotValues(j);
    end
    
    C = [C newCol]; %#ok<AGROW> 
    d = d - C(:,k).^2./C(idx, k);  % Update diagonal.
    [mx, idx] = max( d );   % maximum stays on diagonal because pos. def.
    
end
sz = size(C, 2);       % numerical rank of H.
C = C * spdiags( sqrt( pivotValues ), 0, sz, sz );    % share out scaling.

% Upper-triangular Toeplitz matrix in A = D1(T.*H)D2:
T_row = Lambda2( [1 1:N-1] );
T_row(1) = gamma(alpha-gam+1)/(alpha-gam);
Z = [T_row(1) ; zeros(N-2, 1)];

% Fast Toeplitz matrix multiply. This is the optimized since this is the
% majority of the cost of the code:
b = ifft( bsxfun(@times, fft( bsxfun(@times, C, D2*v), 2*N-1, 1 ),...
    fft( [Z ; 0 ; T_row(end:-1:2)'] ) ) );
c_jac = D1 * ( C .* b(1:N,:) * ones(sz,1) );

% Fix the first entry of the output.
Matrow1 = gamma(gam+beta+2)./gamma(beta+1).*diag(D2)'.*T_row.*vals(1:N)';
c_jac(1) = Matrow1 * v;
c_jac(1) = c_jac(1) + v(1); 

end

function v = UpJacobiConversion(v, a, b)
%UPJACOBICONVERSION   Convert Jacobi (alpha,beta) -> (alpha,beta+1).
%
% UPJACOBICONVERSION(V,A,B) converts coefficients in a Jacobi P^(A,B)
% basis to coefficients in a Jacobi P^(A,B+1) basis in O(length(V)) operations.

% The conversion is upper-triangular and bidiagonal:
[N, M] = size(v);
d1 = [1 (a+b+2)/(a+b+3) (a+b+3:a+b+N)./(a+b+5:2:a+b+2*N-1)].';  % diagonal
d2 = (((a+1:a+N-1)./(a+b+3:2:a+b+2*N-1))).';   % superdiagonal
D1 = spdiags(d1,0,N,N); D2 = spdiags(d2,0,N-1,N-1);
v = D1*v + [ D2*v(2:end,:); zeros(1,M) ];      % Apply conversion matrix
end

function v = DownJacobiConversion(v, a, b)
%DOWNJACOBICONVERSION   Convert Jacobi (alpha,beta+1) -> (alpha,beta).
%
% DOWNJACOBICONVERSION(V, A, B) converts coefficients in a Jacobi P^(A,B+1)
% basis to coefficients in a Jacobi P^(A,B) basis in O(lengt(V)) operations.

% Inversion of UPJACOBICONVERSION, invert upper-triangular and bidiagonal
% matrix fast.

N = length(v);

% first row of inverse up-conversion
topRow = [1 (a+1)/(a+b+2) (a+1)/(a+b+2)*cumprod((a+2:a+N-1)./(a+b+3:a+b+N))];
topRow = (-1).^(0:N-1).*topRow;

% %  Compute S\u in O(N) operations % %
% Sum up topRow.*u' in a numerically stable way.
vecsum = fliplr(cumsum(fliplr(topRow.*v.')));

% Apply inverse up-conversion to u.
ratios = ((a+b+5:2:a+b+2*N-1)./(a+b+3:a+b+N)).*(1./topRow(3:end));
ratios = [ 1 -(a+b+3)/(a+1) ratios];
v = ratios.*vecsum;
v = v.';
end

function v = RightJacobiConversion(v, a, b)
%RIGHTJACOBICONVERSION Convert Jacobi (alpha,beta) -> (alpha+1,beta).
%
% RIGHTJACOBICONVERSION(V,A,B) converts coefficients in a Jacobi P^(A,B)
% basis to coefficients in a Jacobi P^(A+1,B) basis in O(length(V))
% operations.

% Use the reflection formula:
v(2:2:end) = -v(2:2:end);
v = UpJacobiConversion(v, b, a);
v(2:2:end) = -v(2:2:end);
end

function v = LeftJacobiConversion(v, a, b)
%LEFTJACOBICONVERSION Convert Jacobi (alpha+1,beta) -> (alpha,beta).
%
% LEFTJACOBICONVERSION(N,A,B) converts coefficients in a Jacobi P^(A+1,B)
% basis to coefficients in a Jacobi P^(A,B) basis in O(N) operations.

v(2:2:end) = -v(2:2:end);
v = DownJacobiConversion(v, b, a);   % Use the reflection formula.
v(2:2:end) = -v(2:2:end);
end