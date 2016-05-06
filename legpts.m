function [x, w, v, t] = legpts(n, int, meth)
%LEGPTS    Legendre points and Gauss-Legendre Quadrature Weights.
%   LEGPTS(N) returns N Legendre points X in (-1,1).
%
%   [X, W] = LEGPTS(N) returns also a row vector W of weights for Gauss-Legendre
%   quadrature.
%
%   [X, W] = LEGPTS(N, INTERVAL) scales the nodes and weights for the finite
%   interval INTERVAL.
%
%   [X, W, V] = LEGPTS(N) or [X, W, V] = LEGPTS(N, D) returns additionally a
%   column vector V of weights in the barycentric formula corresponding to the
%   points X. The weights are scaled so that max(abs(V)) = 1.
%
%   LEGPTS(N, INTERVAL, METHOD) or LEGPTS(N, METHOD) allows the user to select
%   which method to use.
%    METHOD = 'REC' uses the recurrence relation for the Legendre polynomials
%     and their derivatives to perform Newton iteration on the WKB approximation
%     to the roots. Default for N < 100.
%    METHOD = 'ASY' uses Bogaert's fast algorithm based upon asymptotic 
%     formulae, which is fast and accurate for large N. Default for N >= 100.
%    METHOD = 'GW' uses the traditional Golub-Welsch eigenvalue method,
%     which is maintained mostly for historical reasons.
%
%   [X, W, V, T] = LEGPTS(...) returns also the arccos of the nodes (scaled to
%   lie in [-1, 1] if the INTERVAL argument is used), T = acos(X).  In some
%   situations (in particular with 'ASY') these can be computed to a much
%   better relative precision than X.
%
% See also CHEBPTS, JACPTS, LOBPTS, RADAUPTS, HERMPTS, LAGPTS, and TRIGPTS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTES AND REFERENCES:
%  'GW' by Nick Trefethen, March 2009 - algorithm adapted from [1].
%  'REC' by Nick Hale, July 2011.
%  'ASY' algorithm by Bogaert [2]. Matlab code by Nick Hale, July 2014.
%
%  References:
%   [1] G. H. Golub and J. A. Welsch, "Calculation of Gauss quadrature rules",
%       Math. Comp. 23, 221-230, 1969.
%   [2] I. Bogaert, "Iteration-free computation of Gauss-Legendre quadrature
%       nodes and weights", SIAM J. Sci. Comput., 36(3), A1008-A1026, 2014.
%   [3] A. Glaser, X. Liu and V. Rokhlin, "A fast algorithm for the calculation 
%       of the roots of special functions", SIAM J. Sci. Comput., 2007.
%   [4] N. Hale and A. Townsend, "Fast computation of Gauss-Jacobi quadrature 
%       nodes and weights", SIAM J. Sci. Comput., 2012.
%
%  Historical note:
%   March 2009 - GW [1] algorithm.
%   April 2009 - GLR [3] added for N >= 129.
%     Feb 2011 - REC for N < 129, GLR for large N.
%     Aug 2012 - HT [4] replaces GLR for N >= 129.
%    July 2014 - Bogaert's algorithm [2] replaces HT for N >= 100.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defaults
interval = [-1, 1];
method = 'default';
method_set = nargin == 3;


% Check the inputs:
if ( nargin > 1 )
    if ( nargin == 3 )
        % Calling sequence = LEGPTS(N, INTERVAL, METHOD)
        interval = int;
        method = meth;
    elseif ( nargin == 2 )
        if ( ischar(int) )
            % Calling sequence = LEGPTS(N, METHOD)
            method = int;
            method_set = true;
        else
            % Calling sequence = LEGPTS(N, INTERVAL)
            interval = int;
        end
    end
    validStrings = {'default', 'GW', 'ASY', 'REC'};
    if ( ~any(strcmpi(method, validStrings)) )
        if ( strcmpi(method, 'GLR') )
            error('CHEBFUN:legpts:glr', ...
                'The GLR algorithm is no longer supported.');
        end
        error('CHEBFUN:legpts:inputs', ['Unrecognised input string: ', method]);
    end
    if ( numel(interval) > 2 )
        warning('CHEBFUN:legpts:domain',...
            'Piecewise intervals are not supported and will be ignored.');
        interval = interval([1, end]);
    end
end
if ( any(isinf(interval)) )
    error('CHEBFUN:legpts:interval', 'Unbounded intervals are not supported.');
end

% Deal with trivial cases:
if ( n < 0 )
    error('CHEBFUN:legpts:nNegative', ...
        'First input should be a positive number.');
elseif ( n == 0 )   % Return empty vectors if n == 0:
    x = [];
    w = [];
    v = [];
    t = [];
    return
elseif ( n == 1 )
    x = mean(interval);
    w = diff(interval);
    v = 1;
    t = pi/2;
    return
elseif ( n == 2 )
    x0 = [-1 ; 1]/sqrt(3);
    x = diff(interval)/2 * (x0+1) + interval(1); % map from [-1,1] to interval. 
    w = [1 1]*diff(interval)/2;
    v = [1 ; -1];
    t = acos(x0);
    return
end

if ( n <= 20 )
    % Force REC for n <= 20:
    method = 'rec'; % 
    method_set = 1; 
end

% Choose the method:
t = [];
if ( (n < 100 && ~method_set) || strcmpi(method, 'rec') )
    [x, w, v] = rec(n);        % REC (Standard recurrence relation)
elseif ( strcmpi(method, 'GW') )
    [x, w, v] = gw(n);         % GW see [1]
else
    [x, w, v, t] = asy(n, nargout); % ASY see [2]
end

% Normalise the barycentric weights:
if ( nargout > 2 )
    v = abs(v);
    v = v./max(v);
    v(2:2:end) = -v(2:2:end);
end

% Compute a T is one is asked for:
if ( nargout == 4 && isempty(t) )
    t = acos(x);
end

% Rescale to arbitrary finite interval:
if ( ~all(interval == [-1 1]) )
    dab = diff(interval);
    x = (x+1)/2*dab + interval(1);
    w = dab*w/2;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------------ Routines for GW algorithm ---------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, w, v] = gw(n)
beta = .5./sqrt(1-(2*(1:n-1)).^(-2)); % 3-term recurrence coeffs
T = diag(beta, 1) + diag(beta, -1);   % Jacobi matrix
[V,D] = eig(T);                       % Eigenvalue decomposition
x = diag(D);                          % Legendre points
[x, i] = sort(x);                     % Sort
w = 2*V(1,i).^2;                      % Quadrature weights
v = sqrt(1-x.^2).*abs(V(1,i))';       % Barycentric weights

% Enforce symmetry:
ii = 1:floor(n/2);
x = x(ii);
w = w(ii);
vmid = v(floor(n/2) + 1);
v = v(ii);
if ( mod(n, 2) )
    % Odd number.
    x = [x ; 0 ; -x(end:-1:1)];
    w = [w,  2 - sum(2*w), w(end:-1:1)];
    v = [v ; vmid ; v(end:-1:1)];
else
    % Evem number.
    x = [x ; -x(end:-1:1)];
    w = [w, w(end:-1:1)];
    v = [v ; v(end:-1:1)];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----------------------- Routines for REC algorithm ---------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, w, v] = rec(n)

% Asymptotic formula (Tricomi) - only for positive x.
if ( mod(n,2) )
    s = 1;
else
    s = 0;
end
k = ((n+s)/2:-1:1).';
theta = pi*(4*k-1)/(4*n+2);
x = ( 1 - (n-1)/(8*n^3) - 1/(384*n^4)*(39-28./sin(theta).^2) ).*cos(theta);

% Initialise:
Pm2 = 1;
Pm1 = x;
PPm2 = 0;
PPm1 = 1;
dx = inf;
counter = 0;

% Loop until convergence:
while ( norm(dx, inf) > eps && counter < 10 )
    counter = counter + 1;
    for k = 1:n-1,
        P = ((2*k+1)*Pm1.*x-k*Pm2)/(k+1);
        Pm2 = Pm1;
        Pm1 = P;
        PP = ((2*k+1)*(Pm2+x.*PPm1)-k*PPm2)/(k+1);
        PPm2 = PPm1;
        PPm1 = PP;
    end
    % Newton step:
    dx = -P./PP;
    % Newton update:
    x = x + dx;
    % Reinitialise:
    Pm2 = 1;
    Pm1 = x;
    PPm2 = 0;
    PPm1 = 1;
end

% Once more for derivatives:
for k = 1:n-1,
    P = ( (2*k+1)*Pm1.*x - k*Pm2 ) / (k+1);
    Pm2 = Pm1;
    Pm1 = P;
    PP = ( (2*k+1)*(Pm2+x.*PPm1) - k*PPm2 ) / (k+1);
    PPm2 = PPm1;
    PPm1 = PP;
end

% [TOD0]: This relation might prove useful?
%     PP = -n*(x.*P-Pm2)./(1-x.^2);

% Reflect for negative values:
x = [-x(end:-1:1+s) ; x];
ders = [PP(end:-1:1+s) ; PP];

% Quadrature weights:
w = 2./((1-x.^2).*ders.^2)';

% Barycentric weights
v = 1./ders;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------- Routines for ASY algorithm ------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, w, v, t] = asy(n, nout)

% Compute roots of BesselJ(0, x);
m = ceil(n/2);
jk = besselroots(0,m);

% Useful values:
vn = 1./(n + .5);
a = jk*vn;
u = cot(a);
ua = u.*a;
u2 = u.^2;
a2 = a.^2;

% Initialise for storage (so as to only compute once):
Jk2 = [];
u3 = []; a3 = [];
u4 = []; a4 = [];
u5 = []; a5 = [];
u6 = []; a6 = [];

% Nodes:
[x, t] = legpts_nodes();

% Quadrature weights:
if ( nout > 1 )
    w = legpts_weights();
else
    w = [];
end

% Barycentric weights:
if ( nout > 2 )
    % TODO: Compute bary weights via the given formula in [(6.5)-(6.8), 2]?
%     v = legpts_baryweights();
    v = sin(t)./sqrt(2./w);
    v = v./v(end);
else
    v = [];
end

% Use symmetry:
if ( ~mod(n, 2) )
    x = [-x ; x(end:-1:1)];
    w = [w ; w(end:-1:1)].';
    v = [v ; v(end:-1:1)];
    t = [pi-t ; t(end:-1:1)];
else
    x = [-x(1:end-1) ; 0 ; x(end-1:-1:1)];
    w = [w(1:end) ; w(end-1:-1:1)].';
    v = [v(1:end) ; v(end-1:-1:1)];
    t = [pi-t ; t(end-1:-1:1)];
end

    function [x, t] = legpts_nodes()
        % TODO: Include higher-order terms (i.e., F_4 and F_5).
        
        % Assemble coefficients:
        F0 = a; 
        F1 = 1/8*(u.*a-1)./a;
        if ( n < 1e4 )
            a3 = a.^3;
            F2 = 1/384*( 6*a2.*(1+u2) + 25 - u.*(31*u2+33).*a3 )./a3;
        else
            F2 = 0;
        end
        if ( n < 1e3 )
            u4 = u.^4;
            a5 = a.^5;
            R30 = u.*(2595 + 6350*u2 + 3779*u4)/15360;
            R31 = -(31*u2 + 11)/1024; 
            R32 = u/512; 
            R33 = -25/3072; 
            R35 = -1073/5120;
            F3 = R30 + R35./a5 + (1+u2).*(R31./a + R32./a2 + R33./a3);
        else
            F3 = 0;
        end       
        
        % The asymptotic expansion for the roots:
        t = F0 + F1*vn^2 + F2*vn^4 + F3*vn^6;
        
        % Convert to physical space:
        x = cos(t);
    end

    function w = legpts_weights()
        % TODO: Include higher-order terms (i.e., W_4 and W_5).
        
        % Assemble coefficients:
        W0 = 1; 
        W1 = 1/8*(ua + a2 - 1)./a.^2;
        if ( n < 1e4 )
            a3 = a.^3;
            a4 = a2.^2;
            u4 = u.^4;
            W2 = 1/384*( 81 - 31*ua - 3*(1-2*u2).*a2 + 6*u.*a3 - ...
                (27 + 84*u2 + 56*u4).*a4 )./a4;
        else
            W2 = 0;
        end
        if ( n < 1e3 )
            u3 = u.^3;
            u5 = u.^5;
            u6 = u3.^2;
            a5 = a.^5;
            a6 = a3.^2;
            Q30 = 187/96*u4 + 295/256*u2 + 151/160*u6 + 153/1024;
            Q31 = -119/768*u.^3 -35/384*u5 - 65/1024*u;
            Q32 = 5/512 + 7/384*u4 + 15/512*u2; 
            Q33 = u3/512 - 13/1536*u;
            Q34 = -7/384*u.^2 + 53/3072; 
            Q35 = 3749/15360*u; 
            Q36 = -1125/1024;
            W3 = Q30 + Q31./a + Q32./a2 + Q33./a3 + Q34./a4 + ...
                Q35./a5 + Q36./a6;
        else
            W3 = 0;
        end
        
        % Compute the values of Bessel1(j0k)^2:
        Jk2 = bessel12atj0k(m);
        
        % The asymptotic expansion for the weights:
        w = 2./((Jk2/vn.^2).*(a./sin(a)).*(W0 + W1*vn^2 + W2*vn^4 + W3*vn^6));
        
    end
%     function v = legpts_baryweights()
%         % TODO: Include higher-order terms (i.e., V_4 and V_5).
% 
%         % Assemble coefficients:
%         L0 = 1; 
%         L1 = 1/16*(3*(ua).^2 - 3*ua - (a2-1).*(u2+1) - u2 )./a2;
%         if ( n < 1e4 )
%             L2 = 1/512*( 44*ua + 4*ua.^3 + 22*u.*a3 - 4*ua.^4 + 7*ua.^2 - ...
%                 4*u2.*a4 - 8*a2 + 21*a4 - 51 )./a4;
%         else
%             L2 = 0;
%         end
%         if ( n < 1e3 )
%             V30 = -3353/6144*u4 - 671/8192 - 1663/4096*u2 - 3329/15360*u6;
%             V31 = -5/2048*u5 - 47/8192*u - 7/2048*u3;
%             V32 = 1/4096*u4 + 5/8192*u2 + 43/8192;
%             V33 = -227/12288*u - 85/24576*u3; 
%             V34 = -149/8192*u2 + 145/8192; 
%             V35 = -11861/40960*u; 
%             V36 = 4343/8192;
%             L3 = V30 + V31./a + V32./a2 + V33./a3 + V34./a4 + ...
%                 V35./a5 + V36./a6;
%         else
%             L3 = 0;
%         end
% 
%         % The asymptotic expansion for the weights:
%         v = vn.*sqrt(sin(a).^3./(a.*Jk2)).*(L0 + L1*vn^2 + L2*vn^4 + L3*vn^6);
%     end

end

function Jk2 = bessel12atj0k(m)
%BESSEL12ATJ0k   Evaluate besselj(1,x).^2 at roots of besselj(0,x).
% BESSEL12ATJ0K(M) return besselj(1, bessel0Roots(m)).^2.

% Initialise storage:
Jk2 = zeros(m, 1);

% First 10 values are precomputed (using Wolfram Alpha):
Jk2(1:10) = [   0.2695141239419169
                0.1157801385822037
                0.07368635113640822
                0.05403757319811628
                0.04266142901724309
                0.03524210349099610
                0.03002107010305467
                0.02614739149530809
                0.02315912182469139
                0.02078382912226786];
if ( m <= 10 )
    Jk2 = Jk2(1:m);
    return
end

% Use Taylor series of (NIST, 10.17.3) and McMahon's expansion (NIST, 10.21.19):
k = (11:m).';
ak = pi*(k-.25);
ak2inv = (1./ak).^2;
c = [-171497088497/15206400, 461797/1152, -172913/8064, 151/80, -7/24, 0, 2];
% Jk2(k) = 1./(pi*ak).*polyval(c, ak2inv);
Jk2(k) = 1./(pi*ak).*(c(7) + ak2inv.^2.*(c(5) + ak2inv.*(c(4) + ...
    ak2inv.*(c(3) + ak2inv.*(c(2)+ak2inv.*c(1))))));
end
