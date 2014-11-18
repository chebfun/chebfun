
function [x, w, v] = jacpts(n, a, b, int, meth)
%JACPTS  Gauss-Jacobi Quadrature Nodes and Weights.
%   X = JACPTS(N, ALPHA, BETA) returns the N roots of the degree N Jacobi
%   polynomial with parameters ALPHA and BETA (which must both be greater than
%   or equal -1) where the Jacobi weight function is defined by w(x) =
%   (1-x)^ALPHA*(1+x)^BETA.
%
%   [X, W] = JACPTS(N, ALPHA, BETA) returns also a row vector W of weights for
%   Gauss-Jacobi quadrature.
%
%   [X, W, V] = JACPTS(N, ALPHA, BETA) returns additionally a column vector V of
%   weights in the barycentric formula corresponding to the points X.
%
%   JACPTS(N, ALPHA, BETA, INTERVAL, METHOD) or LEGPTS(N, ALPHA, BETA, METHOD)
%   allows the user to select which method to use.
%    METHOD = 'REC' uses the recurrence relation for the Jacobi polynomials
%     and their derivatives to perform Newton iteration on the WKB approximation
%     to the roots. Default for N < 100.
%    METHOD = 'ASY' uses the Hale-Townsend fast algorithm based upon asymptotic
%     formulae, which is fast and accurate. Default for N >= 100.
%    METHOD = 'GW' will use the traditional Golub-Welsch eigenvalue method,
%       which is maintained mostly for historical reasons.
%
%   [X, W, V] = JACPTS(N, ALPHA, BETA, [A, B]) scales the nodes and weights for
%       the finite interval [A,B].
%
%   The cases ALPHA = BETA = -.5 and ALPHA = BETA = .5 correspond to
%   Gauss-Chebyshev nodes and quadrature, and are treated specially (as a closed
%   form of the nodes and weights is available). ALPHA = BETA = 0 calls LEGPTS,
%   which is a more efficient code.
% 
%   When MAX(ALPHA, BETA) > 5 the results may not be accurate. 
%
% See also CHEBPTS, LEGPTS, LOBPTS, RADAUPTS, HERMPTS, LAGPTS, and TRIGPTS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'GW' by Nick Trefethen, March 2009 - algorithm adapted from [1].
% 'REC' by Nick Hale, July 2011
% 'ASY' by Nick Hale & Alex Townsend, May 2012 - see [2].
%
% NOTE: The subroutines DO NOT SCALE the weights (with the exception of GW).
% This is done in the main code to avoid duplication.
%
% References:
%   [1] G. H. Golub and J. A. Welsch, "Calculation of Gauss quadrature
%       rules", Math. Comp. 23:221-230, 1969.
%   [2] N. Hale and A. Townsend, "Fast computation of Gauss-Jacobi 
%       quadrature nodes and weights", SISC, 2012.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defaults:
interval = [-1, 1];
method = 'default';
method_set = 0;

if ( a <= -1 || b <= -1 )
    error('CHEBFUN:jacpts:sizeAB', 'Alpha and beta must be greater than -1')
elseif ( max(a, b) > 5 )
    warning('CHEBFUN:jacpts:largeAB',...
        'MAX(ALPHA, BETA) > 5. Results may not be accurate')
end


% Check inputs:
if ( nargin > 3 )
    if ( nargin == 5 )
        % Calling sequence = LEGPTS(N, INTERVAL, METHOD)
        interval = int;
        method = meth;
        method_set = 1;
    elseif ( nargin == 4 )
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
            error('CHEBFUN:jacpts:glr', ...
                'The GLR algorithm is no longer supported.');
        end
        error('CHEBFUN:jacpts:inputs', ['Unrecognised input string: ', method]);
    end
end

if ( any(isinf(interval)) )  % Inf intervals not yet supported.
    % TODO: How do we scale the weights?
    error('CHEBFUN:jacpts:infinterval', ... 
          'JACPTS() does not yet support infinite intervals');
elseif ( numel(interval) > 2 )
    warning('CHEBFUN:legpts:domain',...
        'Piecewise intervals are not supported and will be ignored.');
    interval = interval([1, end]);
end

% Deal with trivial cases:
if ( n < 0 )
    error('CHEBFUN:jacpts:n', 'First input should be a positive number.');
elseif ( n == 0 )   % Return empty vectors if n == 0:
    x = []; 
    w = []; 
    v = []; 
    return
elseif ( n == 1 )
    x0 = (b-a)/(a+b+2);
    x = diff(interval)/2 * (x0+1) + interval(1); % map from [-1,1] to interval. 
    w = 2^(a+b+1)*beta(a+1, b+1) * diff(interval)/2;
    v = 1;
    return
end

% Special cases:
if ( ~a && ~b )                 % Legendre: alpha = beta = 0
    [x, w, v] = legpts(n, method);
    [x, w] = rescale(x, w, interval, a, b);
    return
elseif ( a == -.5 && b == -.5 ) % Gauss-Chebyshev: alpha = beta = -.5
    [x, ignored, v] = chebpts(n, interval, 1);
    w = repmat(pi/n,1,n);
    [ignored, w] = rescale(x, w, interval, a, b);
    return
elseif ( a == .5 && b == .5 )   % Gauss-Chebyshev2: alpha = beta = .5
    x = chebpts(n+2, 2);     
    x = x(2:n+1);
    w = pi/(n+1)*(1-x.^2).';  
    v = (1-x.^2);  
    v(2:2:end) = -v(2:2:end); 
    [x, w] = rescale(x,w,interval,a,b);
    return
end

% Choose an algorithm:
if ( n < 20 || (n < 100 && ~method_set) || strcmpi(method, 'rec') )
    [x, w, v] = rec(n, a, b); % REC (Recurrence relation)
    
elseif ( strcmpi(method, 'GW') )
    [x, w, v] = gw(n, a, b);  % GW  see [1]
    w = w/sum(w);
    
else
    [x, w, v] = asy(n, a, b); % HT  see [2]
    
end

% Compute the constant for weights:
if ( a && b )
    if ( n > 50 ) % Use asymptotic approximation:
        M = min(20, n-1); 
        C = 1; 
        phi = -a*b/n;
        for m = 1:M
            C = C + phi;
            phi = -(m+a)*(m+b)/(m+1)/(n-m)*phi;
            if ( abs(phi/C) < eps/100 )
                break
            end
        end
        if ( strcmpi(method, 'rec') )
            C = gamma(2+a) * gamma(2+b) / (gamma(2+a+b)*(a+1)*(b+1));
            w = w / sum(w);
        end
    else
        C = gamma(n+a+1)*gamma(n+b+1)/gamma(n+a+b+1)/factorial(n);
    end
    C = 2^(a+b+1)*C;
else
    C = 2^(a+b+1);
end
if ( ~strcmpi(method,'GW') )
    w = C*w; 
end

% Scale the nodes and quadrature weights:
[x, w] = rescale(x, w, interval, a, b);

% Scale the barycentric weights:
v = abs(v); 
v(2:2:end) = -v(2:2:end);
v = v./max(abs(v)); 

end

function [x, w] = rescale(x, w, interval, a, b)
%RESCALE   Rescale nodes and weights to an arbitrary finite interval.
    if ( ~all(interval == [-1, 1]) )
        c1 = .5*sum(interval); 
        c2 = .5*diff(interval);
        w = c2^(a+b+1)*w;
        x = c1 + c2*x;    
    end
end

%% ------------------------- Routines for GW ----------------------------
    
function [x, w, v] = gw(n, a, b)
    ab = a + b;
    ii = (2:n-1)';
    abi = 2*ii + ab;
    aa = [(b - a)/(2 + ab)
          (b^2 - a^2)./((abi - 2).*abi)
          (b^2 - a^2)./((2*n - 2+ab).*(2*n+ab))];
    bb = [2*sqrt( (1 + a)*(1 + b)/(ab + 3))/(ab + 2) ; 
          2*sqrt(ii.*(ii + a).*(ii + b).*(ii + ab)./(abi.^2 - 1))./abi];
    TT = diag(bb,1) + diag(aa) + diag(bb,-1); % Jacobi matrix.
    [V, x] = eig( TT );                       % Eigenvalue decomposition.
    x = diag(x);                              % Jacobi points.
    % Quadrature weights:
    w = V(1,:).^2*( 2^(ab+1)*gamma(a+1)*gamma(b+1)/gamma(2+ab) ); 
    v = sqrt(1-x.^2).*abs(V(1,:))';           % Barycentric weights.
end

%% ------------------------- Routines for REC ---------------------------


function [x, w, v] = rec(n, a, b)
%REC   Compute nodes and weights using recurrrence relation.

   [x1, ders1] = rec_main(n, a, b, 1); % Nodes and P_n'(x)
   [x2, ders2] = rec_main(n, b, a, 0); % Nodes and P_n'(x)
   x = [-x2(end:-1:1) ; x1];
   ders = [ders2(end:-1:1) ; ders1];
   w = 1./((1-x.^2).*ders.^2)';        % Quadrature weights
   v = 1./ders;                        % Barycentric weights
   
end

function [x, PP] = rec_main(n, a, b, flag)
%REC_MAIN   Jacobi polynomial recurrence relation.

% Asymptotic formula (WKB) - only positive x.
if ( flag )
    r = ceil(n/2):-1:1;
else
    r = floor(n/2):-1:1;  
end
C = (2*r+a-.5)*pi/(2*n+a+b+1);
T = C + 1/(2*n+a+b+1)^2 * ((.25-a^2)*cot(.5*C) - (.25-b^2)*tan(.5*C));
x = cos(T).';

% Initialise:
dx = inf; 
l = 0;
% Loop until convergence:
while ( (norm(dx,inf) > sqrt(eps)/1000) && (l < 10) )
    l = l + 1;
    [P, PP] = eval_Jac(x, n, a, b);
    dx = -P./PP; 
    x = x + dx;
end
% Once more for derivatives:
[ignored, PP] = eval_Jac(x, n, a, b);

end

function [P, Pp] = eval_Jac(x, n, a, b)
%EVALJAC   Evaluate Jacobi polynomial and derivative via recurrence relation.

% Initialise:
ab = a + b;
P = .5*(a-b+(ab+2)*x);  
Pm1 = 1; 
Pp = .5*(ab+2);         
Ppm1 = 0; 

% n = 0 case:
if ( n == 0 )
    P = Pm1; 
    Pp = Ppm1; 
end

for k = 1:n-1
    % Useful values:
    A = 2*(k + 1)*(k + ab + 1)*(2*k + ab);
    B = (2*k + ab + 1)*(a^2 - b^2);
    C = prod(2*k + ab + (0:2)');
    D = 2*(k + a)*(k + b)*(2*k + ab + 2);

    % Recurrence:
    Pa1 = ( (B+C*x).*P - D*Pm1 ) / A;
    Ppa1 = ( (B+C*x).*Pp + C*P - D*Ppm1 ) / A;

    % Update:
    Pm1 = P; 
    P = Pa1;  
    Ppm1 =  Pp; 
    Pp = Ppa1;
end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------- Routines for ASY algorithm ------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, w, v] = asy(n, a, b)
%ASY   Compute nodes and weights using asymptotic formulae.

    if ( n <= 20 ) % Use only boundary formula:
        [xbdy, wbdy, vbdy] = asy2(n, a, b, ceil(n/2));  
        [xbdy2, wbdy2, vbdy2] = asy2(n, b, a, floor(n/2));  
        x = [-xbdy2(end:-1:1) ; xbdy];
        w = [wbdy2(end:-1:1), wbdy];
        v = [vbdy2(end:-1:1) ; vbdy];
        return
    end

    % Determine switch between interior and boundary regions:
    nbdy = 10;
    bdyidx1 = n-(nbdy-1):n; 
    bdyidx2 = nbdy:-1:1;

    % Interior formula:
    [x, w, v] = asy1(n,a,b,nbdy);   

    % Boundary formula (right):
    [xbdy, wbdy, vbdy] = asy2(n, a, b, nbdy);  
    x(bdyidx1) = xbdy;  
    w(bdyidx1) = wbdy; 
    v(bdyidx1) = vbdy;
    
    % Boundary formula (left):
    if ( a ~= b )
        [xbdy, wbdy, vbdy] = asy2(n, b, a, nbdy);  
    end
    x(bdyidx2) = -xbdy; 
    w(bdyidx2) = wbdy; 
    v(bdyidx2) = vbdy;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             ASY1 (Interior)                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, w, v] = asy1(n, a, b, nbdy)
% Algorithm for computing nodes and weights in the interior.

    % Approximate roots via asymptotic formula: (Gatteschi and Pittaluga, 1985)
    K = (2*(n:-1:1)+a-.5)*pi/(2*n+a+b+1);
    tt = K + 1/(2*n+a+b+1)^2*((.25-a^2)*cot(.5*K)-(.25-b^2)*tan(.5*K));

    % First half (x > 0):
    t = tt(tt <= pi/2);
    mint = t(end-nbdy+1);
    idx = 1:max(find(t < mint,1)-1, 1);

    dt = inf; j = 0;
    % Newton iteration
    while ( norm(dt,inf) > sqrt(eps)/100 && j < 10 )
        [vals, ders] = feval_asy1(n, a, b, t, idx, 0);  % Evaluate
        dt = vals./ders;                                % Newton update
        t = t + dt;                                     % Next iterate
        j = j + 1;
        dt = dt(idx);
    end
    [vals, ders] = feval_asy1(n, a, b, t, idx, 1);      % Once more for luck
    t = t + vals./ders;                                 % Newton update.

    % Store:
    x = cos(t);
    w = 1./ders.^2;
    v = (sin(t)./ders);

    % Second half (x < 0):
    tmp = a; 
    a = b; 
    b = tmp;
    t = pi - tt(1:(n-length(x)));
    mint = t(nbdy);
    idx = max(find(t > mint, 1), 1):numel(t);

    dt = inf; j = 0;
    % Newton iteration
    while ( norm(dt,inf) > sqrt(eps)/100 && j < 10 )
        [vals, ders] = feval_asy1(n, a, b, t, idx, 0);  % Evaluate.
        dt = vals./ders;                                % Newton update.
        t = t + dt;                                     % Next iterate.
        j = j + 1;
        dt = dt(idx);
    end
    [vals, ders] = feval_asy1(n, a, b, t, idx, 1);      % Once more for luck.
    t = t + vals./ders;                                 % Newton update.

    % Store:
    x = [-cos(t) x].';
    w = [1./ders.^2 w];
    v = [sin(t)./ders v].';

end

% -------------------------------------------------------------------------

function [vals, ders] = feval_asy1(n, a, b, t, idx, flag)
% Evaluate the interior asymptotic formula at x = cos(t).
    
    % Number of terms in the expansion:
    M = 20;

    % Some often used vectors/matrices:
    onesT = ones(1,length(t));
    onesM = ones(M,1);
    MM = transpose(0:M-1);

    % The sine and cosine terms:
    alpha = (.5*(2*n+a+b+1+MM))*onesT .* (onesM*t) - .5*(a+.5)*pi;
    cosA = cos(alpha);
    sinA = sin(alpha);

    if ( flag )
        if ( idx(1) == 1 )
            k = numel(t):-1:1;
        else
            k = 1:numel(t);
        end
        ta = double(single(t));   
        tb = t - ta;
        hi = n*ta;                
        lo = n*tb + (a+b+1)*.5*t;
        pia = double(single(pi));
        pib = -8.742278000372485e-08; %pib = pi - pia;
        dh = (hi - (k-.25)*pia) + lo - .5*a*pia - (k-.25+.5*a)*pib;
        tmp = 0;
        sgn = 1; fact = 1; DH = dh; dh2 = dh.*dh;
        for j = 0:20
            dc = sgn*DH/fact;
            tmp = tmp + dc;
            sgn = -sgn;
            fact = fact*(2*j+3)*(2*j+2);
            DH = DH.*dh2;
            if ( norm(dc,inf) < eps/2000 )
                break
            end
        end
        tmp(2:2:end) = -tmp(2:2:end);
        tmp = sign(cosA(1,2)*tmp(2))*tmp;
        cosA(1,:) = tmp;
    end

    sinT = onesM*sin(t);
    cosT = onesM*cos(t);
    cosA2 = cosA.*cosT + sinA.*sinT;
    sinA2 = sinA.*cosT - cosA.*sinT;

    one = ones(1,length(t));
    sinT = [one ; cumprod(onesM(2:end)*(.5*csc(.5*t)))];
    cosT = .5*sec(.5*t);

    j = 0:M-2;
    vec = (.5+a+j).*(.5-a+j)./(j+1)./(2*n+a+b+j+2);
    P1 = [1  cumprod(vec)];
    P1(3:4:end) = -P1(3:4:end);
    P1(4:4:end) = -P1(4:4:end);
    P2 = eye(M);
    for l = 0:M-1
        j = 0:(M-l-2);
        vec = (.5+b+j).*(.5-b+j)./(j+1)./(2*n+a+b+j+l+2);
        P2(l+1+(1:length(j)),l+1) = cumprod(vec);
    end
    PHI = repmat(P1,M,1).*P2;

    j = 0:M-2;
    vec = (.5+a+j).*(.5-a+j)./(j+1)./(2*(n-1)+a+b+j+2);
    P1 = [1  cumprod(vec)];
    P1(3:4:end) = -P1(3:4:end);
    P1(4:4:end) = -P1(4:4:end);
    P2 = eye(M);
    for l = 0:M-1
        j = 0:(M-l-2);
        vec = (.5+b+j).*(.5-b+j)./(j+1)./(2*(n-1)+a+b+j+l+2);
        P2(l+1+(1:length(j)),l+1) = cumprod(vec);
    end
    PHI2 = repmat(P1,M,1).*P2;

    S = 0; S2 = 0;
    SC = sinT;
    for m = 0:M-1

        l = 0:2:m;
        phi = PHI(m+1,l+1);
        dS1 = phi*SC(l+1,:).*cosA(m+1,:);

        phi2 = PHI2(m+1,l+1);
        dS12 = phi2*SC(l+1,:).*cosA2(m+1,:);

        l = 1:2:m;
        phi = PHI(m+1,l+1);
        dS2 = phi*SC(l+1,:).*sinA(m+1,:);

        phi2 = PHI2(m+1,l+1);
        dS22 = phi2*SC(l+1,:).*sinA2(m+1,:);

        if m > 10 && norm(dS1(idx) + dS2(idx),inf) < eps/100, break, end

        S = S + dS1 + dS2;
        S2 = S2 + dS12 + dS22;

        SC(1:m+1,:) = bsxfun(@times,SC(1:m+1,:),cosT);
    end

    % Constant out the front:
    dsa = .5*(a^2)/n; dsb = .5*(b^2)/n; dsab = .25*(a+b)^2/n;
    ds = dsa + dsb - dsab; s = ds; j = 1; 
    dsold = ds; % to fix a = -b bug.
    while ( (abs(ds/s) + dsold) > eps/10 )
        dsold = abs(ds/s);
        j = j+1;
        tmp = -(j-1)/(j+1)/n;
        dsa = tmp*dsa*a;
        dsb = tmp*dsb*b;
        dsab = .5*tmp*dsab*(a+b);
        ds = dsa + dsb - dsab;
        s = s + ds;
    end
    p2 = exp(s)*sqrt(2*pi)*sqrt((n+a)*(n+b)/(2*n+a+b))/(2*n+a+b+1);
    g = [1 1/12 1/288 -139/51840 -571/2488320 163879/209018880 ...
         5246819/75246796800 -534703531/902961561600 ...
         -4483131259/86684309913600 432261921612371/514904800886784000];
    f = @(z) sum(g.*[1 cumprod(ones(1,9)./z)]);
    C = p2*(f(n+a)*f(n+b)/f(2*n+a+b))*2/pi;
    C2 = C*(a+b+2*n).*(a+b+1+2*n)./(4*(a+n).*(b+n));

    vals = C*S;
    S2 = C2*S2;

    % Use relation for derivative:
    ders = (n*(a-b-(2*n+a+b)*cos(t)).*vals + 2*(n+a)*(n+b)*S2)/(2*n+a+b)./sin(t);
    denom = 1./real(sin(t/2).^(a+.5).*cos(t/2).^(b+.5));
    vals = vals.*denom;
    ders = ders.*denom;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             ASY2 (Boundary)                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, w, v] = asy2(n, a, b, npts)
% Algorithm for computing nodes and weights near the boundary.

% Use Newton iterations to find the first few Bessel roots:
smallK = min(30, npts);
jk = besselRoots(a, min(npts, smallK));
% Use asy formula for larger ones (See NIST 10.21.19, Olver 1974 p247)
if ( npts > smallK )
    mu = 4*a^2;
    a8 = 8*((length(jk)+1:npts).'+.5*a-.25)*pi;
    jk2 = .125*a8-(mu-1)./a8 - 4*(mu-1)*(7*mu-31)/3./a8.^3 - ...
          32*(mu-1)*(83*mu.^2-983*mu+3779)/15./a8.^5 - ...
          64*(mu-1)*(6949*mu^3-153855*mu^2+1585743*mu-6277237)/105./a8.^7;
    jk = [jk ; jk2];
end
jk = real(jk(1:npts));

% Approximate roots via asymptotic formula: (see Olver 1974)
phik = jk/(n + .5*(a + b + 1));
t = phik + ((a^2-.25)*(1-phik.*cot(phik))./(8*phik) - ...
    .25*(a^2-b^2)*tan(.5*phik))/(n + .5*(a + b + 1))^2;

% Only first half (x > 0):
if ( any(t > 1.1*pi/2) )
    warning('CHEBFUN:jacpts:asy2:theta', ...
        'Theta > pi/2. Result may be inaccurate.'); 
end

% Compute higher terms:
[tB1, A2, tB2, A3] = asy2_higherterms(a, b, t, n);

dt = inf; j = 0;
% Newton iteration:
while ( norm(dt,inf) > sqrt(eps)/200 && j < 10)
    [vals, ders] = feval_asy2(n, t, 1); % Evaluate via asymptotic formula.
    dt = vals./ders;                    % Newton update.
    t = t + dt;                         % Next iterate.
    j = j + 1;
end
[vals, ders] = feval_asy2(n, t, 1);     % Evaluate via asymptotic formula.
dt = vals./ders;                        % Newton update
t = t + dt;    
    
% flip:
t = t(npts:-1:1); 
ders = ders(npts:-1:1);
% vals = vals(npts:-1:1);

% Revert to x-space:
x = cos(t);      
w = (1./ders.^2).';   
v = sin(t)./ders;

    function [vals, ders] = feval_asy2(n, t, flag)
    % Evaluate the boundary asymptotic formula at x = cos(t).
    
        % Useful constants:
        rho = n + .5*(a + b + 1); 
        rho2 = n + .5*(a + b - 1);
        A = (.25 - a^2);       
        B = (.25 - b^2);
        
        % Evaluate the Bessel functions:
        Ja = besselj(a, rho*t, 0);
        Jb = besselj(a + 1, rho*t, 0);
        Jbb = besselj(a + 1, rho2*t, 0);
        if ( ~flag )
            Jab = besselj(a, rho2*t, 0);
        else
            % In the final step, perform accurate evaluation
            Jab = besselTaylor(-t, rho*t, a);
        end

        % Evaluate functions for recurrsive definition of coefficients:
        gt = A*(cot(t/2)-(2./t)) - B*tan(t/2);
        gtdx = A*(2./t.^2-.5*csc(t/2).^2) - .5*B*sec(t/2).^2;
        tB0 = .25*gt;
        A10 = a*(A+3*B)/24;
        A1 = gtdx/8 - (1+2*a)/8*gt./t - gt.^2/32 - A10;
        % Higher terms:
        tB1t = tB1(t); 
        A2t = A2(t); 

        % VALS:
        vals = Ja + Jb.*tB0/rho + Ja.*A1/rho^2 + Jb.*tB1t/rho^3 + Ja.*A2t/rho^4;
        % DERS:
        vals2 = Jab + Jbb.*tB0/rho2 + Jab.*A1/rho2^2 + Jbb.*tB1t/rho2^3 + Jab.*A2t/rho2^4;
        
        % Higher terms (not needed for n > 1000).
        tB2t = tB2(t); A3t = A3(t);
        vals = vals + Jb.*tB2t/rho^5 + Ja.*A3t/rho^6;
        vals2 = vals2 + Jbb.*tB2t/rho2^5 + Jab.*A3t/rho2^6;
        
        % Constant out the front (Computed accurately!)
        ds = .5*(a^2)/n;
        s = ds; jj = 1;
        while abs(ds/s) > eps/10
            jj = jj+1;
            ds = -(jj-1)/(jj+1)/n*(ds*a);
            s = s + ds;
        end
        p2 = exp(s)*sqrt((n+a)/n)*(n/rho)^a;
        g = [1 1/12 1/288 -139/51840 -571/2488320 163879/209018880 ...
             5246819/75246796800 -534703531/902961561600 ...
             -4483131259/86684309913600 432261921612371/514904800886784000];
        f = @(z) sum(g.*[1 cumprod(ones(1,9)./z)]);
        C = p2*(f(n+a)/f(n))/sqrt(2);

        % Scaling:
        valstmp = C*vals;
        denom = sin(t/2).^(a+.5).*cos(t/2).^(b+.5);
        vals = sqrt(t).*valstmp./denom;

        % Relation for derivative:
        C2 = C*n/(n+a)*(rho/rho2)^a;
        ders = (n*(a-b-(2*n+a+b)*cos(t)).*valstmp + 2*(n+a)*(n+b)*C2*vals2)/(2*n+a+b);
        ders = ders.*(sqrt(t)./(denom.*sin(t)));
        
    end

end

function jk = besselRoots(nu, m)
% BESSELROOTS(NU, M) returns the first M roots of besselj(nu, x).
    
% Find an approximation:
jk = zeros(m,1); 
m1 = 3;
if ( nu == 0 )
    xs = 2.404825557695773;
elseif ( nu > 0 )
    % See Hethcote 1970:
    xs = nu + 1.8557*nu^(1/3);
else
    nu1 = nu + 1;
    % See Piessens 1984:
    xs = 2*sqrt(nu+1)*(1 + nu1/4 - 7*nu1^2/96 + 49*nu1^3/1152 - 8363*nu1/276480);
    m1 = min(max(2*ceil(abs(log10(nu1))), 3), m);
end

% The first root:
jk(1) = besselNewton(nu, xs);
if ( m == 1 )
    return
end
% The second root:
jk(2) = besselNewton(nu, jk(1)+.9*pi);
if ( m == 2) 
    return
end
% Some more roots:
for k = 3:m1
    jk(k) = besselNewton(nu, jk(k-1)+.99*pi);
end
% The rest
for k = m1+1:m
    jk(k) = besselNewton(nu, jk(k-1)+pi);
end

end

function jk = besselNewton(nu, jk)
% BESSELNEWTON(NU, JK)   Find roots of Bessel function using Newton iteration.

dx = inf; j = 0;
while ( dx > sqrt(eps)/1000 && j < 20 )
    u = besselj(nu, jk, 0);
    du = besselj(nu-1, jk, 0) - nu/jk*u;
    dx = u./du;
    jk = jk - dx;
    j = j + 1;
end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [TODO]: The codes below are only here temporarily.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ja = besselTaylor(t, z, a)
%BESSELTAYLOR    Accurate evaluation of Bessel function J_A for asy2. (See [2].)
% BESSELTAYLOR(T, Z, A) evaluates J_A(Z+T) by a Taylor series expansion about Z. 

npts = numel(t);
kmax = min(ceil(abs(log(eps)/log(norm(t, inf)))), 30);
H = bsxfun(@power, t, 0:kmax).';
% Compute coeffs in Taylor expansions about z (See NIST 10.6.7)
[nu, JK] = meshgrid(-kmax:kmax, z);
Bjk = besselj(a + nu, JK, 0);
nck = abs(pascal(floor(1.25*kmax), 1)); nck(1,:) = []; % nchoosek
AA = [Bjk(:,kmax+1), zeros(npts, kmax)];
fact = 1;
for k = 1:kmax
    sgn = 1;
    for l = 0:k
        AA(:,k+1) = AA(:,k+1) + sgn*nck(k,l+1)*Bjk(:,kmax+2*l-k+1);
        sgn = -sgn;
    end
    fact = k*fact;
    AA(:,k+1) = AA(:,k+1)/2^k/fact;
end
% Evaluate Taylor series:
Ja = zeros(npts, 1);
for k = 1:npts
    Ja(k,1) = AA(k,:)*H(:,k);
end
end

function [tB1, A2, tB2, A3, tB3, A4] = asy2_higherterms(alph, bet, theta, n)
%ASY2_HIGHERTERMS   Higher-order terms for boundary asymptotic series.
% Compute the higher order terms in asy2 boundary formula. See [2]. 

% These constants are more useful than alph and bet:
A = (0.25 - alph^2);
B = (0.25 - bet^2);

% For now, just work on half of the domain:
c = max(max(theta), .5);
if ( n < 30 )
    N = ceil(40 - n);
elseif ( n >= 30 && c > pi/2-.5)
    N = 15;
else
    N = 10;
end
Nm1 = N - 1;

% Scaled 2nd-kind Chebyshev points and barycentric weights:
t = .5*c*( sin(pi*(-Nm1:2:Nm1)/(2*Nm1)).' + 1 );
v = [.5 ; ones(Nm1,1)];
v(2:2:end) = -1;
v(end) = .5*v(end);

% The g's:
g = A*(cot(t/2)  -2./t) - B*tan(t/2);
gp = A*(2./t.^2 - .5*csc(t/2).^2) - .5*(.25-bet^2)*sec(t/2).^2;
gpp = A*(-4./t.^3 + .25*sin(t).*csc(t/2).^4) - 4*B*sin(t/2).^4.*csc(t).^3;
g(1) = 0; gp(1) = -A/6-.5*B; gpp(1) = 0;

% B0:
B0 = .25*g./t;
B0p = .25*(gp./t-g./t.^2);
B0(1) = .25*(-A/6-.5*B);
B0p(1) = 0;

% A1:
A10 = alph*(A+3*B)/24;
A1 = .125*gp - (1+2*alph)/2*B0 - g.^2/32 - A10;
A1p = .125*gpp - (1+2*alph)/2*B0p - gp.*g/16;
A1p_t = A1p./t;
A1p_t(1) = -A/720 - A^2/576 - A*B/96 - B^2/64 - B/48 + alph*(A/720 + B/48);

% Make f accurately: (Taylor series approx for small t)
fcos = B./(2*cos(t/2)).^2;
f = -A*(1/12 + t.^2/240+t.^4/6048 + t.^6/172800 + t.^8/5322240 + ...
    691*t.^10/118879488000 + t.^12/5748019200 + ...
    3617*t.^14/711374856192000 + 43867*t.^16/300534953951232000);
idx = t > .5;
ti = t(idx);
f(idx) = A.*(1./ti.^2 - 1./(2*sin(ti/2)).^2);
f = f - fcos;

% Integrals for B1: (Note that N isn't large, so we don't need to be fancy).
C = chebcolloc2.cumsummat(N)*(.5*c);
D = chebcolloc2.diffmat(N)*(2/c);
I = (C*A1p_t);
J = (C*(f.*A1));

% B1:
tB1 = -.5*A1p - (.5+alph)*I + .5*J;
tB1(1) = 0;
B1 = tB1./t;
B1(1) = A/720 + A^2/576 + A*B/96 + B^2/64 + B/48 + ...
    alph*(A^2/576 + B^2/64 + A*B/96) - alph^2*(A/720 + B/48);

% A2:
K = C*(f.*tB1);
A2 = .5*(D*tB1) - (.5+alph)*B1 - .5*K;
A2 = A2 - A2(1);

if ( nargout < 3 )
    % Make function for output
    tB1 = @(theta) bary(theta,tB1,t,v);
    A2 = @(theta) bary(theta,A2,t,v);
    return
end

% A2p:
A2p = D*A2;
A2p = A2p - A2p(1);
A2p_t = A2p./t;
% Extrapolate point at t = 0:
w = pi/2-t(2:end);
w(2:2:end) = -w(2:2:end);
w(end) = .5*w(end);
A2p_t(1) = sum(w.*A2p_t(2:end))/sum(w);

% B2:
tB2 = -.5*A2p - (.5+alph)*(C*A2p_t) + .5*C*(f.*A2);
B2 = tB2./t;
% Extrapolate point at t = 0:
B2(1) = sum(w.*B2(2:end))/sum(w);

% A3:
K = C*(f.*tB2);
A3 = .5*(D*tB2) - (.5+alph)*B2 - .5*K;
A3 = A3 - A3(1);

if ( nargout < 6 )
    % Make function for output
    tB1 = @(theta) bary(theta, tB1, t, v);
    A2 = @(theta) bary(theta, A2, t, v);
    tB2 = @(theta) bary(theta, tB2, t, v);
    A3 = @(theta) bary(theta, A3, t, v);
    return
end

% A2p:
A3p = D*A3;
A3p = A3p - A3p(1);
A3p_t = A3p./t;
% Extrapolate point at t = 0:
w = pi/2-t(2:end);
w(2:2:end) = -w(2:2:end);
A3p_t(1) = sum(w.*A3p_t(2:end))/sum(w);

% B2:
tB3 = -.5*A3p - (.5+alph)*(C*A3p_t) + .5*C*(f.*A3);
B3 = tB3./t;
% Extrapolate point at t = 0
B3(1) = sum(w.*B3(2:end))/sum(w);

% A3:
K = C*(f.*tB3);
A4 = .5*(D*tB3) - (.5+alph)*B3 - .5*K;
A4 = A4 - A4(1);

% Make function for output:
tB1 = @(theta) bary(theta, tB1, t, v);
A2 = @(theta) bary(theta, A2, t, v);
tB2 = @(theta) bary(theta, tB2, t, v);
A3 = @(theta) bary(theta, A3, t, v);
tB3 = @(theta) bary(theta, tB3, t, v);
A4 = @(theta) bary(theta, A4, t, v);

end



