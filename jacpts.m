function [x,w,v] = jacpts(n,a,b,varargin)
%JACPTS  Gauss-Jacobi Quadrature Nodes and Weights.
%  X = JACPTS(N,ALPHA,BETA) returns the N roots of the degree N Jacobi 
%       polynomial with parameters ALPHA and BETA (which must both be 
%       greater than or equal -1) where the Jacobi weight function is 
%       defined by w(x) = (1-x)^ALPHA*(1+x)^BETA.
%
%  [X,W] = JACPTS(N,ALPHA,BETA) returns also a row vector W of weights for 
%       Gauss-Jacobi quadrature.
%
%  [X,W,V] = JACPTS(N,ALPHA,BETA) returns additionally a column vector V of 
%       weights in the barycentric formula corresponding to the points X.
%
%  [X,W,V] = JACPTS(N,ALPHA,BETA,METHOD) allows choice in which method is 
%       used. METHOD = 'GW' will use the traditional Golub-Welsch
%       eigenvalue method, which is best suited for when N is small. METHOD
%       = 'GLR' will use the Glaser-Liu-Rokhlin fast algorithm, which is
%       much faster for large N. By default JACPTS will use 'GW' when N <
%       128.
%
%  [X,W,V] = JACPTS(N,ALPHA,BETA,[A,B]) scales the nodes and weights for the 
%       finite interval [A,B].
%
%  The cases ALPHA = BETA = -.5 and ALPHA = BETA = .5 correspond to
%  Gauss-Chebyshev nodes and quadrature, and are treated specially (as a
%  closed form of the nodes and weights is available). ALPHA = BETA = 0
%  calls LEGPTS, which is a more efficient code.
%
%  See also legpts and chebpts.

%  Copyright 2011 by The University of Oxford and The Chebfun Developers. 
%  See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

%  'GW' by Nick Trefethen, March 2009 - algorithm adapted from [1].
%  'GLR' by Nick Hale, April 2009 - algorithm adapted from [2].
%  'REC' by Nick Hale, July 2011
%  'ASY' by Nick Hale & Alex Townsend, May 2012 - see [3].
%
%   NOTE: The subroutines DO NOT SCALE the weights (with the exception of
%   GW). This is done in the main code to avoid duplication.
%
%  References:
%   [1] G. H. Golub and J. A. Welsch, "Calculation of Gauss quadrature
%       rules", Math. Comp. 23:221-230, 1969, 
%   [2] A. Glaser, X. Liu and V. Rokhlin, "A fast algorithm for the 
%       calculation of the roots of special functions", SIAM Journal  
%       on Scientific Computing", 29(4):1420-1438:, 2007.
%   [3] N. Hale and A. Townsend, "Fast computation of Gauss-Jacobi 
%       quadrature nodes and weights",In preparation, 2012.


% Defaults
interval = [-1,1];
method = 'default';
method_set = 0;

if n < 0
    error('CHEBFUN:jacpts:n','First input should be a positive number.');
end
if a <= -1 || b <= -1,
    error('CHEBFUN:jacpts:SizeAB','alpha and beta must be greater than -1')
end
if nargout > 1 && any(isinf(interval)) % infinite intervals not yet supported
    error('CHEBFUN:jacpts:infinterval', ... % (How do we scale the weights?) 
    'jacpts does not yet support infinite intervals');
end

% Check inputs
if nargin > 3
    if isa(varargin{1},'double') && length(varargin{1}) == 2
        interval = varargin{1};
    elseif isa(varargin{1},'domain')
        interval = varargin{1}.ends;
    elseif isa(varargin{1},'char')
        method = varargin{1}; 
        method_set = 1;
    end
    if length(varargin) == 2,
        if isa(varargin{2},'double') && length(varargin{2}) == 2
            interval = varargin{2};
        elseif isa(varargin{1},'domain')
            interval = varargin{2}.ends;
        elseif isa(varargin{2},'char')
            method = varargin{2}; 
            method_set = 1;
        end
    end
end

% Deal with trivial cases
if n < 0
    error('CHEBFUN:legpts:n', 'First input should be a positive number.');
elseif n == 0   % Return empty vector if n == 0
    x = []; w = []; v = []; return
elseif n == 1 || n == 2
    method = 'rec'; method_set = 1; % force REC for n = 1 or 2
end

% Special cases
if ~(a || b)                % Legendre: alpha = beta = 0
    [x, w, v] = legpts(n,varargin{:});
    [x, w] = rescale(x,w,interval,a,b);
    return
elseif a == -.5 && b == -.5 % Gauss-Chebyshev: alpha = beta = -.5
    [x, ignored, v] = chebpts(n,interval,1); %#ok<*ASGLU>
    w = repmat(pi/n,1,n);
    [x, w] = rescale(x,w,interval,a,b);
    return
elseif a == .5 && b == .5   % Gauss-Chebyshev2: alpha = beta = .5
    x = chebpts(n+2,2);     x = x(2:n+1);
    w = pi/(n+1)*(1-x.^2);  w = w';
    if nargout == 3, v = (1-x.^2);  v(2:2:end) = -v(2:2:end); end
    [x, w] = rescale(x,w,interval,a,b);
    return
end

% Choose an algorithm
if (n < 100 && ~method_set) || any(strcmpi(method,{'fastsmall','rec'}))
    [x, w, v] = rec(n,a,b); % REC ('fastsmall' is for backward compatibiilty)
elseif strcmpi(method,'GW')
    [x, w, v] = gw(n,a,b);  % GW  see [1]
elseif strcmpi(method,'GLR')  
    [x, w, v] = glr(n,a,b); % GLR see [2]
    w = w/sum(w);
else
    [x, w, v] = asy(n,a,b); % HT  see [3]
end

% Constant for weights
if a && b
    if n > 50
        M = min(20,n-1); C = 1; phi = -a*b/n;
        for m = 1:M
            C = C + phi;
            phi = -(m+a)*(m+b)/(m+1)/(n-m)*phi;
            if abs(phi/C) < eps/100, break, end
        end
        if any(strcmpi(method,{'rec','fastsmall','GLR'}))
            C = gamma(2+a)*gamma(2+b)/(gamma(2+a+b)*(a+1)*(b+1));
            w = w/sum(w);
        end
    else
        C = gamma(n+a+1)*gamma(n+b+1)/gamma(n+a+b+1)/factorial(n);
    end
    C = 2^(a+b+1)*C;
else
    C = 2^(a+b+1);
end
if ~strcmpi(method,'GW'), w = C*w; end

% Scaling
[x w] = rescale(x,w,interval,a,b);
v = abs(v); v(2:2:end) = -v(2:2:end);
v = v./max(abs(v)); 

end

function [x w] = rescale(x,w,interval,a,b)
    % rescale to arbitrary interval
    if all(interval == [-1 1])
        % Nothing to do
        return
    elseif ~any(isinf(interval))
        % finite interval
        c1 = .5*sum(interval); 
        c2 = .5*diff(interval);
        w = c2^(a+b+1)*w;
        x = c1+c2*x;        
    else
        % infinite interval (not yet supported)
        m = maps(fun,{'unbounded'},interval); % use default map
        if nargout > 1
            w = w.*m.der(x.');
        end
        x = m.for(x);
        x([1 end]) = interval([1 end]);
    end
end

%% ------------------------- Routines for GW ----------------------------
    
function [x w v] = gw(n,a,b)
    ab = a + b;
    ii = (2:n-1)';
    abi = 2*ii + ab;
    aa = [(b - a)/(2 + ab)
          (b^2 - a^2)./((abi - 2).*abi)
          (b^2 - a^2)./((2*n - 2+ab).*(2*n+ab))];
    bb = [2*sqrt( (1 + a)*(1 + b)/(ab + 3))/(ab + 2) ; 
          2*sqrt(ii.*(ii + a).*(ii + b).*(ii + ab)./(abi.^2 - 1))./abi];
    TT = diag(bb,1) + diag(aa) + diag(bb,-1); % Jacobi matrix
    [V x] = eig( TT );                        % Eigenvalue decomposition
    x = diag(x);                              % Jacobi points
    % Quadrature weights
    w = V(1,:).^2*( 2^(ab+1)*gamma(a+1)*gamma(b+1)/gamma(2+ab) ); 
    v = sqrt(1-x.^2).*abs(V(1,:))';           % Barycentric weights
end

%% ------------------------- Routines for REC ---------------------------


function [x w v] = rec(n,a,b)

   [x1 ders1] = rec_main(n,a,b,1);          % Nodes and P_n'(x)
   [x2 ders2] = rec_main(n,b,a,0);          % Nodes and P_n'(x)
   x = [-x2(end:-1:1) ; x1];
   ders = [ders2(end:-1:1) ; ders1];
   w = 1./((1-x.^2).*ders.^2)';        % Quadrature weights
   v = 1./ders;                        % Barycentric weights
   
end

function [x PP] = rec_main(n,a,b,flag)

    % Asymptotic formula (WKB) - only positive x.
    if flag, r = ceil(n/2):-1:1;
    else     r = floor(n/2):-1:1;  end
    C = (2*r+a-.5)*pi/(2*n+a+b+1);
    T = C + 1/(2*n+a+b+1)^2*((.25-a^2)*cot(.5*C)-(.25-b^2)*tan(.5*C));
    x = cos(T).';

    % Initialise
    dx = inf; l = 0;
    % Loop until convergence
    while norm(dx,inf) > sqrt(eps)/1000 && l < 10
        l = l + 1;
        [P PP] = eval_Jac(x,n,a,b);
        dx = -P./PP; x = x + dx;
    end
    % Once more for derivatives
    [ignored PP] = eval_Jac(x,n,a,b);

end

%% ------------------------- Routines for GLR ---------------------------

function [x w v] = glr(n,a,b)
    % Main routine for GLR. 

    if abs(a)<=.5 && abs(b)<=.5 % use asymptotic formula
        
        % Choose a root near the middle
        r = ceil(n/2); 
        C = (2*r+a-.5)*pi/(2*n+a+b+1);
        T = C + 1/(2*n+a+b+1)^2*((.25-a^2)*cot(.5*C)-(.25-b^2)*tan(.5*C));
        x1 = cos(T);
        % Make accurate using Newton
        [u up] = eval_Jac(x1,n,a,b);
        while abs(u) > eps
            x1 = x1 - u/up;
            [u up] = eval_Jac(x1,n,a,b);      
        end

        [rootsl dersl] = alg1_Jac_as(n,x1,up,a,b,1);     % Get roots to the left
        if a ~= b 
            [rootsr dersr] = alg1_Jac_as(n,x1,up,a,b,0); % To the right
        else
            rootsr = -rootsl; % Use symmetry.
            dersr = dersl;
            if rootsl(1) > 0, rootsl(1) = []; dersl(1) = []; end
        end

        x = [rootsl(end:-1:2) ; rootsr];
        ders = [dersl(end:-1:2) ; dersr];
        
    else % Start at the middle (more accurate than end)

        [x1 up] = alg3_Jac(n,0,a,b);                  % Find a root

        [rootsl dersl] = alg1_Jac(n,x1,up,a,b,1);     % Get roots to the left
        if a ~= b 
            [rootsr dersr] = alg1_Jac(n,x1,up,a,b,0); % To the right
        else
            rootsr = -rootsl(1+mod(n,2):end);         % Use symmetry.
            dersr = dersl(1+mod(n,2):end);
            if rootsl(1) > 0, rootsl(1) = []; dersl(1) = []; end
        end

        x = [rootsl(end:-1:2) ; rootsr];    
        ders = [dersl(end:-1:2) ; dersr];
    end

    w = 1./((1-x.^2).*ders.^2)';        % Quadrature weights
    v = 1./ders;                        % Barycentric weights

end

% ---------------------------------------------------------------------

function [roots ders] = alg1_Jac(n,x1,up,a,b,flag)
    % Workhorse of GLR routine. Step from root to root.

    ab = a + b;
    if flag, sgn = -1; else sgn = 1; end
    N = n-1;
    roots = [x1 ; zeros(N,1)]; ders = [up ; zeros(N,1)];
    m = 30; % number of terms in Taylor expansion
    hh1 = ones(m+1,1); u = zeros(1,m+1); up = zeros(1,m+1);
    for j = 1:N
        x = roots(j);
        h = rk2_Jac(sgn*pi/2,-sgn*pi/2,x,n,a,b) - x;
        if abs(x+h) > 1, roots(j+1:end) = []; ders(j+1:end) = []; return, end
        
        % scaling
        M = 1/h;

        % recurrence relation  (scaled)
        u(1) = 0;   u(2) = ders(j)/M;  up(1) = u(2); up(m+1) = 0;    
        for k = 0:m-2
            u(k+3) = ((x*(2*k+ab+2)+a-b)*u(k+2)/M + ...
                (k*(k+ab+1)-n*(n+ab+1))*u(k+1)/M^2/(k+1))./((1-x.^2)*(k+2));
            up(k+2) = (k+2)*u(k+3)*M;
        end

        % flip for more accuracy in inner product calculation
        u = u(m+1:-1:1);
        up = up(m+1:-1:1);

        % Newton iteration
        hh = hh1; hh(end) = M;    step = inf;  l = 0; 
        while (abs(step) > eps) && (l < 10)
            l = l + 1;
            step = (u*hh)/(up*hh);
            h = h - step;
            hh = [M;cumprod(M*h+zeros(m,1))]; % powers of h (This is the fastest way!)
            hh = hh(end:-1:1);
        end

        if abs(h) < eps, roots(j+1:end) = []; ders(j+1:end) = []; return, end

        % update
        roots(j+1) = x + h;
        ders(j+1) = up*hh;   
    end

end
% -------------------------------------------------------------------------

function [roots ders] = alg1_Jac_as(n,x1,up,a,b,flag) 
    % if |a|<=.5 && |b|<=.5 use asymptotic formula
    ab = a + b;

    % Approximate roots via asymptotic formula
    nx1 = ceil(n/2);
    if flag, r = (nx1+1):n; else r = (nx1-1):-1:1; end
    C = (2*r+a-.5)*pi/(2*n+ab+1);
    T = C + 1/(2*n+a+b+1)^2*((.25-a^2)*cot(.5*C)-(.25-b^2)*tan(.5*C));
    roots = [x1 ; cos(T).']; ders = [up ; zeros(length(T),1)];

    m = 30; % number of terms in Taylor expansion
    hh1 = ones(m+1,1); u = zeros(1,m+1); up = zeros(1,m+1);
    for j = 1:length(r)
        x = roots(j); % previous root

        % initial approx (via asymptotic foruma)
        h = roots(j+1) - x;
        % scaling
        M = 1/h;

        % recurrence relation for Jacobi polynomials (scaled)
        cx = (1-x.^2);
        u(1) = 0;   u(2) = ders(j)/M;  up(1) = u(2); up(m+1) = 0;    
        for k = 0:m-2
            u(k+3) = ((x*(2*k+ab+2)+a-b)*u(k+2)/M + ...
                (k*(k+ab+1)-n*(n+ab+1))*u(k+1)/M^2/(k+1))./(cx*(k+2));
            up(k+2) = (k+2)*u(k+3)*M;
        end

        % flip for more accuracy in inner product calculation
        u = u(m+1:-1:1);  up = up(m+1:-1:1);

        % Newton iteration
        hh = hh1; hh(end) = M;    step = inf;  l = 0; 
        while (abs(step) > eps) && (l < 10)
            l = l + 1;
            step = (u*hh)/(up*hh);
            h = h - step;
            hh = [M;cumprod(M*h+zeros(m,1))]; % powers of h (This is the fastest way!)
            hh = hh(end:-1:1); % flip for more accuracy in inner product calculation
        end

        % update
        roots(j+1) = x + h;
        ders(j+1) = up*hh;    
    end
end

% ---------------------------------------------------------------------

function [x1 d1] = alg3_Jac(n,xs,a,b)
    % A Newton iteration to find the first root

    [u up] = eval_Jac(xs,n,a,b);
    theta = atan((1-xs.^2)*up/sqrt(n*(n+a+b+1)*(1-xs.^2))/u);
    x1 = rk2_Jac(theta,-pi/2,xs,n,a,b);

    for k = 1:10
        [u up] = eval_Jac(x1,n,a,b);
        x1 = x1 - u/up;
    end

    [ignored d1] = eval_Jac(x1,n,a,b);

end

% -------------------------------------------------------------------------

function [P Pp] = eval_Jac(x,n,a,b)
    % Evaluate Jacobi polynomial and derivative via recurrence relation
    
    ab = a + b;
    P = .5*(a-b+(ab+2)*x);  Pm1 = 1; 
    Pp = .5*(ab+2);         Ppm1 = 0; 

    if n == 0; P = Pm1; Pp = Ppm1; end

    for k = 1:n-1
        A = 2*(k+1)*(k+ab+1)*(2*k+ab);
        B = (2*k+ab+1)*(a^2-b^2);
        C = prod(2*k+ab+(0:2)');
        D = 2*(k+a)*(k+b)*(2*k+ab+2);

        Pa1 = ( (B+C*x).*P - D*Pm1 ) / A;
        Ppa1 = ( (B+C*x).*Pp + C*P - D*Ppm1 ) / A;

        Pm1 = P; P = Pa1;  Ppm1 =  Pp; Pp = Ppa1;
    end
    
end

% -------------------------------------------------------------------------

function x = rk2_Jac(t,tn,x,n,a,b)
    % March along ODE solution using Runge-Kutta
    
    ab = a + b;
    m = 10; h = (tn-t)/m;
    for j = 1:m
        f1 = (1-x.^2);
        k1 = -4*h*f1./(4*sqrt(n*(n+ab+1)*f1) + (b-a-(ab+1)*x).*sin(2*t));
        t = t+h;
        f2 = (1-(x+k1).^2);
        k2 = -4*h*f2./(4*sqrt(n*(n+ab+1)*f2) + (b-a-(ab+1)*(x+k1)).*sin(2*t));
        x = x+.5*real(k1+k2);
    end
    
end

%% ------------------------- Routines for ASY ---------------------------

function [x w v] = asy(n,a,b)

    if n <= 20 % use only boundary formula
        [xbdy wbdy vbdy] = asy2(n,a,b,ceil(n/2));  
        [xbdy2 wbdy2 vbdy2] = asy2(n,b,a,floor(n/2));  
        x = [-xbdy2(end:-1:1) ; xbdy];
        w = [wbdy2(end:-1:1)  wbdy];
        v = [vbdy2(end:-1:1) ; vbdy];
        return
    end

    % Determine switch between interior and boundary regions
    nbdy = 10;
    bdyidx1 = n-(nbdy-1):n; 
    bdyidx2 = nbdy:-1:1;

    % Interior
    [x w v] = asy1(n,a,b,nbdy);   

    % Boundary
    [xbdy wbdy vbdy] = asy2(n,a,b,nbdy);  
    x(bdyidx1) = xbdy;  w(bdyidx1) = wbdy; v(bdyidx1) = vbdy;
    if a ~= b
        [xbdy wbdy vbdy] = asy2(n,b,a,nbdy);  
    end
    x(bdyidx2) = -xbdy; w(bdyidx2) = wbdy; v(bdyidx2) = vbdy;
    
end

% -------------------------------------------------------------------------

function [x w v] = asy1(n,a,b,nbdy)

    % Approximate roots via asymptotic formula.
    K = (2*(n:-1:1)+a-.5)*pi/(2*n+a+b+1);
    tt = K + 1/(2*n+a+b+1)^2*((.25-a^2)*cot(.5*K)-(.25-b^2)*tan(.5*K));

    % First half, x > 0
    t = tt(tt<=pi/2);
    mint = t(end-nbdy+1);
    idx = 1:max(find(t<mint,1)-1,1);

    dt = inf; j = 0;
    % Newton iteration
    while norm(dt,inf) > sqrt(eps)/100
        [vals ders] = feval_asy1(n,a,b,t,idx,0); % Evaluate via asymptotic formulae
        dt = vals./ders;                   % Newton update
        t = t + dt;                        % Next iterate
        j = j + 1;
        dt = dt(idx);
        if j > 10, dt = 0; end
    end
    [vals ders] = feval_asy1(n,a,b,t,idx,1); % Once more for luck
    t = t + vals./ders;

    % Store
    x = cos(t);
    w = 1./ders.^2;
    v = (sin(t)./ders);
%     t1 = t;

    % Second half, x < 0
    tmp = a; a = b; b = tmp;
    t = pi - tt(1:(n-length(x)));
    mint = t(nbdy);
    idx = max(find(t>mint,1),1):numel(t);

    dt = inf; j = 0;
    % Newton iteration
    while norm(dt,inf) > sqrt(eps)/100
        [vals ders] = feval_asy1(n,a,b,t,idx,0); % Evaluate via asymptotic formulae
        dt = vals./ders;                  % Newton update
        t = t + dt;                       % Next iterate
        j = j + 1;
        dt = dt(idx);
        if j > 10, dt = 0; end
    end
    [vals ders] = feval_asy1(n,a,b,t,idx,1); % Once more for luck
    t = t + vals./ders;                  % Newton update

    x = [-cos(t) x].';
    w = [1./ders.^2 w];
    v = [sin(t)./ders v].';
%     t = [t t1].';

end

% -------------------------------------------------------------------------

function [vals ders] = feval_asy1(n,a,b,t,idx,flag) %#ok<INUSD>
    
    % Number of terms in the expansion
    M = 20;

    % Some often used vectors/matrices
    onesT = ones(1,length(t));
    onesM = ones(M,1);
    MM = transpose(0:M-1);

    % The sine and cosine terms
    alpha = (.5*(2*n+a+b+1+MM))*onesT .* (onesM*t) - .5*(a+.5)*pi;
    cosA = cos(alpha);
    sinA = sin(alpha);

    % if flag
        if idx(1) == 1
            k = numel(t):-1:1;
        else
            k = 1:numel(t);
        end
        ta = double(single(t));   tb = t - ta;
        hi = n*ta;                lo = n*tb + (a+b+1)*.5*t;
        pia = double(single(pi));
        pib = -8.742278000372485e-08; %pib = pi - pia;
        dh = (hi-(k-.25)*pia)+lo-.5*a*pia-(k-.25+.5*a)*pib;
        tmp = 0;
        sgn = 1; fact = 1; DH = dh; dh2 = dh.*dh;
        for j = 0:20
            dc = sgn*DH/fact;
            tmp = tmp + dc;
            sgn = -sgn;
            fact = fact*(2*j+3)*(2*j+2);
            DH = DH.*dh2;
            if norm(dc,inf) < eps/2000, break, end
        end
        tmp(2:2:end) = -tmp(2:2:end);
        tmp = sign(cosA(1,2)*tmp(2))*tmp;
        cosA(1,:) = tmp;
    % end

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

    % Constant out the front
    dsa = .5*(a^2)/n; dsb = .5*(b^2)/n; dsab = .25*(a+b)^2/n;
    ds = dsa + dsb - dsab; s = ds; j = 1; 
    dsold = ds; % to fix a = -b bug.
    while (abs(ds/s) + dsold) > eps/10 
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

    % Use relation for derivative
    ders = (n*(a-b-(2*n+a+b)*cos(t)).*vals + 2*(n+a)*(n+b)*S2)/(2*n+a+b)./sin(t);
    denom = 1./real(sin(t/2).^(a+.5).*cos(t/2).^(b+.5));
    vals = vals.*denom;
    ders = ders.*denom;

end

% -------------------------------------------------------------------------

function [x w v] = asy2(n,a,b,npts)

smallk = min(30,npts);
% Use GLR for finding the first bessel roots
jk = BesselRoots(a,min(npts,smallk));
% use asy formula for larger ones (See NIST 10.21.19, Olver 74 p247)
if npts > smallk
    mu = 4*a^2;
    a8 = 8*((length(jk)+1:npts).'+.5*a-.25)*pi;
    jk2 = .125*a8-(mu-1)./a8 - 4*(mu-1)*(7*mu-31)/3./a8.^3 - ...
          32*(mu-1)*(83*mu.^2-983*mu+3779)/15./a8.^5 - ...
          64*(mu-1)*(6949*mu^3-153855*mu^2+1585743*mu-6277237)/105./a8.^7;
    jk = [jk ; jk2];
end
jk = real(jk(1:npts));

% Approximate roots via asymptotic formula (see Olver 1974)
phik = jk/(n + .5*(a + b + 1));
t = phik + ((a^2-.25)*(1-phik.*cot(phik))./(8*phik) - ...
    .25*(a^2-b^2)*tan(.5*phik))/(n + .5*(a + b + 1))^2;

% Only first half, x > 0
if any(t > pi/2), warning('jacpts_bdy:theta','Theta > pi/2'); end

% Compute higher terms
[tB1 A2 tB2 A3] = asy2_higherterms(a,b,t,n);

dt = inf; j = 0;
% Newton iteration
while norm(dt,inf) > sqrt(eps)/200
    [vals ders] = feval_asy2(n,t,1);   % Evaluate via asymptotic formula
    dt = vals./ders;                   % Newton update
    t = t + dt;                        % Next iterate
    j = j + 1; if j > 10, dt = 0; end  % Bail
end

[vals ders] = feval_asy2(n,t,1);     % Evaluate via asymptotic formula
dt = vals./ders;                   % Newton update
t = t - dt;    
    
% flip
t = t(npts:-1:1); ders = ders(npts:-1:1);% vals = vals(npts:-1:1);
% Revert to x-space
x = cos(t);      w = (1./ders.^2).';   v = sin(t)./ders;

    function [vals ders] = feval_asy2(n,t,flag)
        
        % Useful constants
        rho = n + .5*(a + b + 1); 
        rho2 = n + .5*(a + b - 1);
        A = (.25-a^2);       B = (.25-b^2);
        
        
        % Evaluate the Bessel functions
        Ja = besselj(a,rho*t,0);
        Jb = besselj(a+1,rho*t,0);
        Jbb = besselj(a+1,rho2*t,0);
        if ~flag
            Jab = besselj(a,rho2*t,0);
        else
            % In the final step, perform accurate evaluation
            Jab = BesselTaylor(-t,rho*t);
        end

        % Evaluate functions for recurrsive definition of coefficients.
        gt = A*(cot(t/2)-(2./t)) - B*tan(t/2);
        gtdx = A*(2./t.^2-.5*csc(t/2).^2) - .5*B*sec(t/2).^2;
        tB0 = .25*gt;
        A10 = a*(A+3*B)/24;
        A1 = gtdx/8 - (1+2*a)/8*gt./t - gt.^2/32 - A10;
        tB1t = tB1(t); A2t = A2(t); % Higher terms

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

        valstmp = C*vals;
        denom = sin(t/2).^(a+.5).*cos(t/2).^(b+.5);
        vals = sqrt(t).*valstmp./denom;

        % Relation for derivative
        C2 = C*n/(n+a)*(rho/rho2)^a;
        ders = (n*(a-b-(2*n+a+b)*cos(t)).*valstmp + 2*(n+a)*(n+b)*C2*vals2)/(2*n+a+b);
        ders = ders.*(sqrt(t)./(denom.*sin(t)));
        
    end

    function Ja = BesselTaylor(t,z)
        kmax = min(ceil(abs(log(eps)/log(norm(t,inf)))),30);
        H = bsxfun(@power,t,0:kmax).';
        % Compute coeffs in Taylor expansions about z (See NIST 10.6.7)
        [nu JK] = meshgrid(-kmax:kmax, z);
        Bjk = besselj(a+nu,JK,0);
        nck = abs(pascal(floor(1.25*kmax),1)); nck(1,:) = []; % nchoosek
        AA = [Bjk(:,kmax+1) zeros(npts,kmax)];
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
        % Evaluate Taylor series
        Ja = zeros(npts,1);
        for k = 1:npts
            Ja(k,1) = AA(k,:)*H(:,k);
        end
    end

end

function jk = BesselRoots(nu,m)
    % Find m roots of besselj(nu,x)
    
    jk = zeros(m,1); foo = 3;
    if nu == 0
        xs = 2.404825557695773;
    elseif nu > 0
        % See Hethcote 1970
        xs = nu + 1.8557*nu^(1/3);
    else
        nu1 = nu + 1;
        % See Piessens 1984
        xs = 2*sqrt(nu+1)*(1 + nu1/4 - 7*nu1^2/96  + 49*nu1^3/1152 - 8363*nu1/276480);
        foo = min(max(2*ceil(abs(log10(nu1))),3),m);
    end

    % The first root
    jk(1) = BesselNewton(nu,xs);
    if m == 1, return, end
    % The second root
    jk(2) = BesselNewton(nu,jk(1)+.9*pi);
    if m == 2, return, end
    % Some more roots
    for k = 3:foo
        jk(k) = BesselNewton(nu,jk(k-1)+.99*pi);
    end
    % The rest
    for k = foo+1:m
        jk(k) = BesselNewton(nu,jk(k-1)+pi);
    end

end

function jk = BesselNewton(nu,jk)
    % Use Newton iterations to find roots

    dx = inf; j = 0;
    while dx > sqrt(eps)/1000;
        u = besselj(nu,jk,0);
        du = besselj(nu-1,jk,0)-nu/jk*u;
        dx = u./du;
        jk = jk - dx;
        j = j+1;
        if j > 20, break; end
    end
    
end

% -------------------------------------------------------------------------

function [tB1 A2 tB2 A3 tB3 A4] = asy2_higherterms(a,b,theta,n)
    % Compute the higher order terms in asy2 boundary formula 

    % The constants a = alpha and b = beta
    A = (.25-a^2); B = (.25-b^2); % These are more useful

    % For now, just work on half of the domain
    % c = pi/2; N = 30;
    c = max(max(theta),.5); 
    if n < 30, N = ceil(20-(n-20)); else N = 10; end
    if n > 30 && c > pi/2-.5, N = 15; end
    N1 = N-1;

    % 2nd-kind Chebyshev points and barycentric weights
    t = .5*c*(sin(pi*(-N1:2:N1)/(2*N1)).'+1);        
    v = [.5 ; ones(N1,1)]; v(2:2:end) = -1; v(end) = .5*v(end);

    % The g's
    g = A*(cot(t/2)-2./t)-B*tan(t/2);
    gp = A*(2./t.^2-.5*csc(t/2).^2)-.5*(.25-b^2)*sec(t/2).^2;
    gpp = A*(-4./t.^3+.25*sin(t).*csc(t/2).^4)-4*B*sin(t/2).^4.*csc(t).^3;
    g(1) = 0; gp(1) = -A/6-.5*B; gpp(1) = 0;

    % B0
    B0 = .25*g./t;
    B0p = .25*(gp./t-g./t.^2);
    B0(1) = .25*(-A/6-.5*B);
    B0p(1) = 0;

    % A1
    A10 = a*(A+3*B)/24;
    A1 = .125*gp - (1+2*a)/2*B0 - g.^2/32 - A10;
    A1p = .125*gpp - (1+2*a)/2*B0p - gp.*g/16;
    A1p_t = A1p./t;
    A1p_t(1) = -A/720-A^2/576-A*B/96-B^2/64-B/48+a*(A/720+B/48);

    % Make f accurately
    fcos = B./(2*cos(t/2)).^2;
    f = -A*(1/12+t.^2/240+t.^4/6048+t.^6/172800+t.^8/5322240 + ...
        691*t.^10/118879488000+t.^12/5748019200+3617*t.^14/711374856192000 + ...
        43867*t.^16/300534953951232000);
    idx = t>.5;
    ti = t(idx);
    f(idx) = A.*(1./ti.^2 - 1./(2*sin(ti/2)).^2);
    f = f - fcos;

    % Integrals for B1
    C = cumsummat(N)*(.5*c); 
    D = diffmat(N)*(2/c);
    I = (C*A1p_t);
    J = (C*(f.*A1));

    % B1
    tB1 = -.5*A1p - (.5+a)*I + .5*J;
    tB1(1) = 0;
    B1 = tB1./t;
    B1(1) = A/720+A^2/576+A*B/96+B^2/64+B/48+a*(A^2/576+B^2/64+A*B/96)-a^2*(A/720+B/48);

    % A2
    K = C*(f.*tB1);
    A2 = .5*(D*tB1) - (.5+a)*B1 - .5*K;
    A2 = A2 - A2(1);

    %[tB1b A2b] = asy2_highertermsc(.1,-.3);

    if nargout < 3
        % Make function for output
        tB1 = @(theta) bary(theta,tB1,t,v);
        A2 = @(theta) bary(theta,A2,t,v);
    end

    % A2p
    A2p = D*A2;
    A2p = A2p - A2p(1);
    A2p_t = A2p./t;
    % Extrapolate point at t = 0
    w = pi/2-t(2:end);
    w(2:2:end) = -w(2:2:end);
    w(end) = .5*w(end);
    A2p_t(1) = sum(w.*A2p_t(2:end))/sum(w);

    % B2
    tB2 = -.5*A2p - (.5+a)*(C*A2p_t) + .5*C*(f.*A2);
    B2 = tB2./t;
    % Extrapolate point at t = 0
    B2(1) = sum(w.*B2(2:end))/sum(w);

    % A3
    K = C*(f.*tB2);
    A3 = .5*(D*tB2) - (.5+a)*B2 - .5*K;
    A3 = A3 - A3(1);

    if nargout < 6
        % Make function for output
        tB1 = @(theta) bary(theta,tB1,t,v);
        A2 = @(theta) bary(theta,A2,t,v);
        tB2 = @(theta) bary(theta,tB2,t,v);
        A3 = @(theta) bary(theta,A3,t,v);
        return
    end

    % A2p
    A3p = D*A3;
    A3p = A3p - A3p(1);
    A3p_t = A3p./t;
    % Extrapolate point at t = 0
    w = pi/2-t(2:end);
    w(2:2:end) = -w(2:2:end);
    A3p_t(1) = sum(w.*A3p_t(2:end))/sum(w);

    % B2
    tB3 = -.5*A3p - (.5+a)*(C*A3p_t) + .5*C*(f.*A3);
    B3 = tB3./t;
    % Extrapolate point at t = 0
    B3(1) = sum(w.*B3(2:end))/sum(w);

    % A3
    K = C*(f.*tB3);
    A4 = .5*(D*tB3) - (.5+a)*B3 - .5*K;
    A4 = A4 - A4(1);

    % Make function for output
    tB1 = @(theta) bary(theta,tB1,t,v);
    A2 = @(theta) bary(theta,A2,t,v);
    tB2 = @(theta) bary(theta,tB2,t,v);
    A3 = @(theta) bary(theta,A3,t,v);
    tB3 = @(theta) bary(theta,tB3,t,v);
    A4 = @(theta) bary(theta,A4,t,v);    

end




