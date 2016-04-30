function pass = test_jac2jac( )
% Test for jac2jac 

rng(0)
v = randn(10,2);
% Run through a load of alpha, beta, gam, and delta: 
alpha = linspace(-.99,1.1,5);
beta = linspace(-.99,1.1,5); 
gam = linspace(-.99,1.1,5);
delta = linspace(-.99,1.1,5); 

tol = 2e-10; 
count = 1; 
ind = [];
% Test accuracy against a direct method: 
for a = alpha 
    for b = beta 
        for g = gam 
            for d = delta 
                exact = cheb2jac_direct( jac2cheb_direct(v, a, b), g, d ); 
                w = jac2jac( v, a, b, g, d );
                err(count) =  norm( exact - w, inf ); 
                ind(count,:) =  [ a, b, g, d ];
                pass(count) = norm( exact - w, inf ) < tol; 
                count = count + 1;
            end
        end
    end
end

% Take a look at err and ind if there is a failure.
if ( all(pass) ) 
    pass = 1; 
else
    pass = 0;
end

% Try something over N = 513: 
N = 513; 
v = randn(N,2);
a = .46; b = -.7; g = .56; d = 1.54; 
exact = cheb2jac_direct( jac2cheb_direct( v, a, b), g, d ); 
w = jac2jac( v, a, b, g, d );
pass(2) = norm( exact - w, inf ) < tol;

% Try alpha + beta = -1: 
N = 513; 
v = randn(N,2);
a = -.6; b = -.4; g = -.65; d = -.45; 
exact = cheb2jac_direct( jac2cheb_direct( v, a, b), g, d ); 
w = jac2jac( v, a, b, g, d );
pass(3) = norm( exact - w, inf ) < 2*tol;

% Try alpha - gamma approx 0: 
N = 513; 
v = randn(N,2);
a = .1; b = -.4; g = .10000001; d = -.4; 
exact = cheb2jac_direct( jac2cheb_direct( v, a, b), g, d ); 
w = jac2jac( v, a, b, g, d );
pass(4) = norm( exact - w, inf ) < tol;

end

function c_jac = cheb2jac_direct(c_cheb, a, b)
%CHEB2LEG_DIRECT   Convert Cheb to Leg coeffs using the 3-term recurrence.
[N, m] = size(c_cheb);              % Number of columns.
N = N - 1;                          % Degree of polynomial.
% Don't let N be too big:
if ( N > 2^11 )
    error('CHEBFUN:cheb2jac:arraySize', ...
        'Maximum transform size (2048) is exceeded.');
end
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

function c_cheb = jac2cheb_direct(c_jac, a, b)
%JAC2CHEB_DIRECT   Convert Leg to Cheb coeffs using the 3-term recurrence.
N = size(c_jac,1) - 1;                 % Degree of polynomial.
% Don't let N be too big:
if ( N > 2^11 )
    error('CHEBFUN:jac2cheb:arraySize', ...
        'Maximum transform size (2048) is exceeded.');
end
if ( N <= 0 ), c_cheb = c_jac; return, end  % Trivial case.
tech = chebtech1;                      % Alternatively use chebtech2.
x = tech.chebpts(N+1);                 % Chebyshev grid (reversed order).
% Make the Jacobi-Chebyshev Vandemonde matrix:
apb = a + b; aa  = a * a; bb  = b * b;
P = zeros(N+1); P(:,1) = 1;    
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
v_cheb = P*c_jac;                      % Values on Chebyshev grid.
c_cheb = tech.vals2coeffs(v_cheb);     % Chebyshev coefficients.
end