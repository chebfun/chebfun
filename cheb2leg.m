function c_leg = cheb2leg( c_cheb, varargin )
%CHEB2LEG   Convert Chebyshev coefficients to Legendre coefficients.
%   C_LEG = CHEB2LEG(C_CHEB) converts the vector C_CHEB of Chebyshev
%   coefficients to a vector C_LEG of Legendre coefficients such that
%       C_CHEB(1)*T0 + ... + C_CHEB(N)*T{N-1} = ...
%           C_LEG(1)*P0 + ... + C_LEG(N)*P{N-1},
%   where P{k} is the degree k Legendre polynomial normalized so that
%   max(|P{k}|) = 1.
%
%   C_LEG = CHEB2LEG(C_CHEB, 'norm') is as above, but with the Legendre
%   polynomials normalized to be orthonormal.
%
%   If C_CHEB is a matrix then the CHEB2LEG operation is applied to each column.
%
% See also LEG2CHEB, CHEB2JAC, JAC2CHEB, JAC2JAC.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%% For developers  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For more information on the algorithm, see Section 5.3 of
%
% A. Townsend, M. Webb, and S. Olver, "Fast polynomial transforms based on
% Toeplitz and Hankel matrices", submitted, 2016.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BRIEF IDEA:
% Let A be the upper-triangular conversion matrix. We observe that A can be
% decomposed as A = D1(T.*H)D2, where D1 and D2 are diagonal, T is
% Toeplitz, and H is a real, symmetric, positive definite Hankel matrix.
% The Hankel part can be approximated, up to an error of tol, by a
% rank O( log N log(1/tol) ) matrix. A low rank approximation is
% constructed via pivoted Cholesky factorization.

[N, n] = size(c_cheb);

normalize = 0;                                  % Default - no normalize.
trans = 0; 
for j = 1:numel(varargin)
    if ( strncmpi(varargin{j}, 'norm', 4) )
        normalize = 1;
    elseif ( strncmpi(varargin{j}, 'trans', 4) )
        trans = 1;
    end
end

if ( trans ) 
    stop 
end

% Trivial case:
if ( N < 2 )
    c_leg = c_cheb;
    return
elseif ( N <= 512 ) 
    % Use direct approach: 
    c_leg =  cheb2leg_direct( c_cheb ); 
    if ( normalize ) 
        c_leg  = bsxfun(@times, c_leg, 1./sqrt((0:(N-1))'+.5) ); 
    end
    return
end

%%%%%%%%%%%%%%%%%%%%%% Initialise  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate the symbol of the Hankel part of M:
% This for-loop is a faster and more accurate way of doing:
% Lambda = @(z) exp(gammaln(z+1/2) - gammaln(z+1));
% vals = Lambda( (0:2*N-1)'/2 );
vals = zeros(2*N,1);
vals(1) = sqrt(pi);
vals(2) = 2/vals(1);
for i = 2:2:2*(N-1)
    vals(i+1) = vals(i-1)*(1-1/i);
    vals(i+2) = vals(i)*(1-1/(i+1));
end

% First row of cheb2leg matrix:
num = 0:N-1;
L_row1 = -num/2.*([ 1 1 vals(1:N-2)' ]./num).*([ 1 vals(1:N-1)' ]./(num+1));
L_row1(2:2:end) = 0; L_row1(1) = 1;

% Calculate SVD of hankel part:
% [ii, jj] = meshgrid( num );
% H =  Lambda( abs(ii + jj - 1)/2 )./(ii+jj+1);
% H = H(2:end,2:end);

%%%%%%%%%%%%%%%%%%  Pivoted Cholesky algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate quasi-SVD of hankel part by using Cholesky factorization.
% Find the numerical rank and pivot locations. This is equivalent to
% Cholesky factorization on the matrix A with diagonal pivoting, except
% here only the diagonal of the matrix is updated.
num = (1:N-1)';
d = vals(2:2:2*N-2).*(num.^2./((2:2:2*N-2)'+1)); % diagonal of Hankel matrix.
pivotValues = [];    % store Cholesky pivots.
C = [];              % store Cholesky columns.
tol = 1e-14*log(N);  % Tolerance of low rank approx.
k = 0;
[mx, idx] = max( d );% max on diagonal, searching for first Cholesky pivot.
while ( mx > tol )
    
    k = k + 1;
    
    pivotValues = [pivotValues ; 1./mx];
    
    newCol = vals(idx+1:idx+N-1).*(num.*num(idx)./(idx + (1:N-1)' + 1));
    
    %C = [  C,  H(:, idx) ];
    for j = 1 : k - 1
        newCol = newCol - C(:, j) * C(idx, j) * pivotValues(j);
    end
    C = [ C,  newCol ];
    %C(:, k) = newCol;
    
    d = d - C(:,k).^2./C(idx, k);  % Update diagonal.
    [mx, idx] = max( d );          % Find next pivot.
    
end
sz = size(C, 2);        % numerical rank of H.
C = C * spdiags( sqrt( pivotValues ), 0, sz, sz );    % share out scaling.

% Second row of Toeplitz part:
T_row2 = [ 0 0 vals(1:N-3)']./(num'-1);
T_row2(1) = 0; T_row2(2:2:end) = 0;
Z = [T_row2(1) ; zeros(N-2, 1)];

% Diagonal part of the matrix:
d = sqrt(pi)/2./vals(3:2:2*N);

%%%%%%%%%%%%%%%%%%  Multiply  D1(T.*H)D2  %%%%%%%%%%%%%%%%%%%%%%%%%%
c_leg = zeros(N, n); 

if ( n == 1 )
    b = ifft( bsxfun(@times, fft( bsxfun(@times, C, -c_cheb(2:end)), 2*N-2 ),...
                                    fft( [Z ; 0 ; T_row2(end:-1:2)'] ) ) );
    c_leg(2:end) = ((num+1/2)./num).*((C.* b(1:N-1,:))*ones(sz,1));
    c_leg(2:end) = c_leg(2:end) + d.*c_cheb(2:end);
    c_leg(1) = L_row1 * c_cheb;
else
    a = fft( [Z ; 0 ; T_row2(end:-1:2)'] );
    for k = 1 : n
    b = ifft( bsxfun(@times, fft( bsxfun(@times, C, -c_cheb(2:end,k)), 2*N-2, 1 ), a), [], 1 );
    c_leg(2:end,k) = ((num+1/2)./num).*((C.* b(1:N-1,:))*ones(sz,1));
    c_leg(2:end,k) = c_leg(2:end,k) + d.*c_cheb(2:end,k);
    c_leg(1,k) = L_row1 * c_cheb(:,k);
    end
end

% Normalize:
if ( normalize ),
    c_leg  = bsxfun(@times, c_leg, 1./sqrt((0:N-1)'+.5) );
end

end


function c_leg = cheb2leg_direct(c_cheb)
%CHEB2LEG_DIRECT   Convert Cheb to Leg coeffs using the 3-term recurrence.
[N, m] = size(c_cheb);              % Number of columns.
N = N - 1;                          % Degree of polynomial.
x = cos(.5*pi*(0:2*N).'/N);         % 2*N+1 Chebyshev grid (reversed order).
f = dct1([c_cheb ; zeros(N, m)]);   % Values on 2*N+1 Chebyshev grid.
w = chebtech2.quadwts(2*N+1).';     % Clenshaw-Curtis quadrature weights.
Pm2 = 1; Pm1 = x;                   % Initialise.
L = zeros(2*N+1, N+1);              % Vandermonde matrix.
L(:,1:2) = [1+0*x, x];              % P_0 and P_1.
for k = 1:N-1 % Recurrence relation:
    P = (2-1/(k+1))*Pm1.*x - (1-1/(k+1))*Pm2;    
    Pm2 = Pm1; Pm1 = P; 
    L(:,2+k) = P;
end
scale = (2*(0:N).'+1)/2;            % Scaling in coefficients [NIST, (18.3.1)].
c_leg = bsxfun(@times, L.'*(bsxfun(@times, f ,w)), scale); % Legendre coeffs.
end

function v = dct1(c)
%DCT1   Compute a (scaled) DCT of type 1 using the FFT. 
% DCT1(C) returns T(X)*C, where X = cos(pi*(0:N)/N) and T(X) = [T_0, T_1, ...,
% T_N](X) (where T_k is the kth 1st-kind Chebyshev polynomial), and N =
% length(C) - 1;

c([1,end],:) = 2*c([1,end],:);              % Scale.
v = chebfun.dct(c, 1);                      % DCT-I.

end