function p = legpoly(n, dom, normalize, method)
%LEGPOLY   Legendre polynomials.
%   P = LEGPOLY(N) computes a CHEBFUN of the Legendre polynomial of degree N on
%   the interval [-1,1]. N can be a vector of integers, in which case the output
%   is an array-valued CHEBFUN.
%
%   P = LEGPOLY(N, D) computes the Legendre polynomials as above, but on the
%   interval given by the domain D, which must be bounded.
%
%   P = LEGPOLY(N, D, 'norm') or P = LEGPOLY(N, 'norm') normalises so that
%   integral(P(:,j).*P(:,k)) = delta_{j,k}.
%
%   For N <= 1000 LEGPOLY uses a weighted QR factorisation of a 2*(N+1) x
%   2*(N+1) Chebyshev Vandermonde matrix. For scalar N > 1000 (or a short
%   vector) it uses the LEG2CHEB method and for a vector of N with any entry >
%   1000 it uses the standard recurrence relation. This default can be
%   overwritten by passing a fourth input LEGPOLY(N, D, NORM, METHOD), where
%   METHOD is 1, 2, or 3 respectively.
%
% See also CHEBPOLY, LEGPTS, and LEG2CHEB.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse input:
methodIsSet = false;
if ( isempty(n) )
    p = chebfun; 
    return
end
if ( nargin < 2 || isempty(dom) )
    dom = [-1, 1];
end
if ( nargin < 3 || isempty(normalize) )
    normalize = 0; 
end
if ( ischar(dom) )
    if ( nargin == 3 )
        method = normalize;
        if ( ~isempty(method) )
            methodIsSet = true;
        end
    end
    normalize = dom;
    dom = [-1,1]; 
end
if ( strncmp(normalize, 'norm', 4) )
    normalize = 1;
elseif ( ischar(normalize) )
    normalize = 0; 
end
if ( (nargin == 4) && ~isempty(method) )
    methodIsSet = true;
end
    
% Unbounded domains aren't supported/defined.
if ( any(isinf(dom)) )
    error('CHEBFUN:legpoly:infdomain', ...
        'Legendre polynomials are not defined over an unbounded domain.');
end

% Force a CHEBTECH basis.
defaultPref = chebfunpref();
pref = defaultPref;
tech = feval(pref.tech);
if ( ~isa(tech, 'chebtech') )
    pref.tech = @chebtech2;
end

% Useful values:
nMax = max(n);
nMax1 = nMax + 1;
domIn = dom;
dom = dom([1 end]);

% Determine which method:
if ( ~methodIsSet && nMax > 1000 )
    % Use LEG2CHEB():
    method = 3;
elseif ( ~methodIsSet )
    % Use QR orthogonalization of Chebyshev polynomials for moderate nmax.
    method = 2;
end

% If the user wants most of the Legendre polynomials then it is faster to use
% the recurrence: TODO: "most" is most than 20% of the [0:nMax]. Should this
% percentage vary with nMax?
if ( ~methodIsSet && (method == 3) && (nMax < 5000) && (numel(n) > nMax/5) ) 
    method = 1; % Do not go above 5000 with the recurrence.
end

switch method
    case 1 % Recurrence
        
        [aa, bb, cc] = unique(n);      %#ok<ASGLU>
        P = zeros(nMax1, length(n));   % Initialise storage
        x = chebpts(nMax1, 2);         % Chebyshev points
        L0 = ones(nMax1, 1); L1 = x;   % P_0 and P_1
        ind = 1;                       % Initialise counter
        for k = 2:nMax+2,              % The recurrence relation (k = degree)
            if ( aa(ind) == k-2 )
                if ( normalize )
                    invnrm = sqrt((2*k-3)/diff(dom));
                    P(:,ind) = L0*invnrm;
                else
                    P(:,ind) = L0;
                end
                ind = ind + 1;
            end
            tmp = L1;
            L1 = (2-1/k)*x.*L1 - (1-1/k)*L0;
            L0 = tmp;
        end
        C = chebtech2.vals2coeffs(P(:,cc));       % Convert to coefficients          
    case 2 % QR

        pts = 2*nMax1;              % Expand on Chebyshev grid of twice the size
        [ignored, w] = chebpts(pts, 2);   % Grab the Clenshaw-Curtis weights
        theta = pi*(pts-1:-1:0)'/(pts-1);
        A = cos(theta*(0:nMax));                  % Vandemonde-type matrix
        D = spdiags(sqrt(w(:)), 0, pts, pts);     % C-C quad weights
        Dinv = spdiags(1./sqrt(w(:)), 0, pts, pts);
        [Q, ignored] = qr(D*A, 0);                % Weighted QR
        P = Dinv*Q;
        if ( normalize )
            PP = P(:,n+1) * diag(sqrt(2/diff(dom)) * sign(P(end,n+1)));
        else
            PP = P(:,n+1) * diag(1./P(end,n+1));
        end
        C = chebtech2.vals2coeffs(PP);            % Convert to coefficients
        C(nMax1+1:end,:) = [];                    % Trim coefficients > nMax
        
    case 3 % LEG2CHEB
        
        c_leg = zeros(nMax+1, numel(n));
        c_leg(n+1,:) = eye(numel(n));             % Legendre coefficients          
        if ( normalize )
            C = leg2cheb(c_leg, 'norm');  % Chebyshev coefficients
        else
            C = leg2cheb(c_leg);          % Chebyshev coefficients
        end
    
end

% Construct CHEBFUN from coeffs:
p = chebfun(C, dom, pref, 'coeffs');              

if ( numel(domIn) > 2 )
    p = restrict(p, domIn);
end

% Adjust orientation:
if ( size(n, 1) > 1 )
   p = p.'; 
end

end
