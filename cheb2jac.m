function c_jac = cheb2jac( c_cheb, alpha, beta ) 
%CHEB2JAC   Convert Chebyshev coefficients to Jacobi coefficients. 
%   C_JAC = CHEB2LEG(C_CHEB, A, B) converts the vector C_CHEB of Chebyshev
%   coefficients to a vector C_JAC of Jacobi coefficients such that 
%    C_CHEB(1)*T_0(x) + ... + C_CHEB(N)*T{N-1}(x) = ...
%           C_LEG(1)*P_0^{(A,B)}(x) + ... + C_LEG(N)*P{N-1}^{(A,B)}(x),
%   where P_k^{(A,B)} is the degree k Jacobi polynomial corresponding to the
%   weight function w(x) = (1-X)^A * (1+X)^B.
%
% See also JAC2CHEB, JAC2JAC.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%% For developers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For more information on the algorithm for N>512, see Section 5.3 of 
% 
% A. Townsend, M. Webb, and S. Olver, "Fast polynomial transforms based on 
% Toeplitz and Hankel matrices", submitted, 2016. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(c_cheb, 1);    % Length of coefficients. 

if ( alpha == 0 && beta == 0 ) 
   % Use cheb2leg for alpha = beta = 0:
    
    c_jac = cheb2leg_new( c_cheb ); 
    
elseif ( alpha == -.5 && beta == -.5 ) 
    % Undo scaling: 
    
    % Convert T_n -> P_n^(-1/2,1/2): 
    scl = [1 cumprod((1/2:1/2+N-2)./(1:N-1))]';   % P_n^(-1/2,-1/2)(1)
    c_jac = spdiags( scl, 0, N, N) \ c_cheb;  
    
elseif ( N <= 512 ) 
    
    c_jac = cheb2jac_direct(c_cheb, alpha, beta);
    
else
    % Convert Chebyshev to Jacobi (-1/2,-1/2) and then call jac2jac: 
    
    % Convert T_n -> P_n^(-1/2,1/2): 
    scl = [1 cumprod((1/2:1/2+N-2)./(1:N-1))]';   % P_n^(-1/2,-1/2)(1)
    c_jac = spdiags( scl, 0, N, N) \ c_cheb; 
    % Now convert P_n^(-1/2,-1/2) -> P_n^(alpha,beta): 
    c_jac = jac2jac( c_jac, -1/2, -1/2, alpha, beta );

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% DIRECT METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c_jac = cheb2jac_direct(c_cheb, a, b)
%CHEB2JAC_DIRECT   Convert Cheb to Jacobi (A,B) coeffs using direct method.

[N, m] = size(c_cheb);              % Number of columns.
N = N - 1;                          % Degree of polynomial.
% Don't let N be too big:

if ( N <= 0 ), c_jac = c_cheb; return, end % Trivial case.
f = chebtech2.coeffs2vals([c_cheb ; zeros(N, m)]); % Values on 2*N+1 Cheb grid.
% 2*N+1 Chebyshev grid (reversed order) and Clenshaw-Curtis-Jacobi weights:
[w, x] = ccjQuadwts(2*N+1, a, b); 

% Make the Jacobi-Chebyshev Vandemonde matrix:
apb = a + b; aa  = a * a; bb  = b * b;
P = zeros(2*N+1, N+1); P(:,1) = 1;    
P(:,2) = 0.5*(2*(a + 1) + (apb + 2)*(x - 1));   
for k = 2:N
    k2 = 2*k;
    k2apb = k2 + apb;
    q1 =  k2*(k + apb)*(k2apb - 2);
    q2 = (k2apb - 1)*(aa - bb);
    q3 = (k2apb - 2)*(k2apb - 1)*k2apb;
    q4 =  2*(k + a - 1)*(k + b - 1)*k2apb;
    P(:,k+1) = ((q2 + q3*x).*P(:,k) - q4*P(:,k-1)) / q1;
end

% Scaling:
% NN = (0:N)';
% scale = 2^(a+b+1)*gamma(NN+a+1).*gamma(NN+b+1) ./ ...
%     ((2*NN+a+b+1).*gamma(NN+a+b+1).*factorial(NN))
scale = zeros(N+1, 1);
scale(1) = beta(a+1, b+1);
for n = 0:N-1
    scale(n+2) = (2*n+a+b+1)*(n+a+1)*(n+b+1) / ...
        ((n+1)*(2*n+a+b+3)*(n+a+b+1))*scale(n+1);
end
scale = 2^(a+b+1)*scale;

% Jacobi coefficients:
c_jac = bsxfun(@times, P.'*(bsxfun(@times, f , w.')), 1./scale); 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% DCT METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [w, x] = ccjQuadwts(n, a, b)
%CCJQUADWTS   Clenshaw-Curtis-Jacobi quadrature weights.
%   [W, X] = CCJQUADWTS(N, A, B) returns the N-point Clenshaw-Curtis-Jacobi
%   quadrature nodes, X = CHEBPTS(N), and weights, W, corresponding to the
%   weight function w(t) = (1-t)^A * (1+t)^B on the interval [-1,1].

% TODO: Move this to somewhere more sensible / accessible.

if ( a == b && a == 0 ) % Clenshaw-Curtis

    c = 2./[1, 1-(2:2:(n-1)).^2];          % Standard Chebyshev moments
    c = [c, c(floor(n/2):-1:2)];           % Mirror for DCT via FFT 
    w = ifft(c);                           % Interior weights
    w([1, n]) = w(1)/2;                    % Boundary weights

elseif ( a == b )       % Gegenbauer
    
    l = a + .5;                            % Gegenbauer parameter
    g0 = gamma(l+.5)*sqrt(pi)/gamma(l+1);
    k = 1:floor((n-1)/2); 
    c = g0*[1, cumprod((k-l-1)./(k+l))];   % Chebyshev moments for (1-x)^a(1+x)^b
    c = [c, c(floor(n/2):-1:2)];           % Mirror for DCT via FFT 
    w = ifft(c);                           % Interior weights
    w([1, n]) = w(1)/2;                    % Boundary weights
    
else                    % Jacobi
    
    c = [1, (a-b)/(a+b+2), zeros(1, n-2)]; % Initialise moments
    for r = 1:n % Recurrence relation for 3F2([r, -r, b +1 ; .5, a+b+2] ; 1 ):
        c(r+2) = - (2*(b-a)*c(r+1) + (a+b+2-r)*c(r)) / (a+b+2+r);
    end
    c = 2^(a+b+1)*gamma(a+1)*gamma(b+1)/gamma(a+b+2) * c; % Moments (with const)
    v = ifft([c(1:n), c(n-1:-1:2)]);       % Mirror for DCT via FFT 
    w = [v(1), 2*v(2:n-1), v(n)];          % Rescale interior weights

end

if ( nargout > 1 )
    x = chebtech2.chebpts(n);              % 2nd-kind Chebyshev points.
end

end
