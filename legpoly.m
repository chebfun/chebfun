function p = legpoly(n, dom, normalize, method)
%LEGPOLY   Legendre polynomials.
%   P = LEGPOLY(N) computes a chebfun of the Legendre polynomial of degree N on
%   the interval [-1,1]. N can be a vector of integers, in which case the output
%   is an array-valued CHEBFUN.
%
%   P = LEGPOLY(N, D) computes the Legendre polynomials as above, but on the
%   interval given by the domain D, which must be bounded.
%
%   P = LEGPOLY(N, D, 'norm') or P = LEGPOLY(N, 'norm') normalises so that
%   integrate(P(:,j).*P(:,k)) = delta_{j,k}.
%
%   For N <= 1000 LEGPOLY uses a weighted QR factorisation of a 2*(N+1) x
%   2*(N+1) Chebyshev Vandermonde matrix. For scalar N > 1000 it uses the
%   LEG2CHEB method and for a vector of N with any entry > 1000 it uses the
%   standard recurrence relation. This default can be overwritten by passing a
%   fourth input LEGPOLY(N, D, NORM, METHOD), where METHOD is 1, 2, or 3
%   respectively.
%
%   NOTE, LEGPOLY() will always return a CHEBFUN whose underlying 'tech' is a
%   CHEBTECH2 object.
%
% See also CHEBPOLY, LEGPTS, and LEG2CHEB.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: This code needs a test.

% Parse input:
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
    normalize = dom;
    dom = [-1,1]; 
end
if ( strncmp(normalize, 'norm', 4) )
    normalize = 1;
elseif ( ischar(normalize) )
    normalize = 0; 
end
    
% Unbounded domains aren't supported/defined.
if ( any(isinf(dom)) )
    error('CHEBFUN:legpoly:infdomain', ...
        'Legendre polynomials are not defined over an unbounded domain.');
end

% Useful values:
nMax = max(n);
nMax1 = nMax + 1;
dom = dom([1 end]); % [TODO]: Add support for breakpoints?

% Determine which method
if ( nargin == 4 )
    % Method has been passed as an input
elseif ( nMax > 1000 )
    if ( numel(n) == 1 )
        % Use LEG2CHEB():
        method = 3;
    else
        % Use three-term recurrence (possible loss of orthogonality)
        method = 1;
    end
else
    % Use QR orthogonalization of Chebyshev polynomials for moderate nmax
    method = 2;
end
if ( ~isscalar(n) && method == 3 )
    % Can only use method 3 for scalar inputs.
    method = 2;
end

switch method
    case 1 % Recurrence
        
        [aa, bb, cc] = unique(n);      %#ok<ASGLU>
        P = zeros(nMax1, length(n));   % Initialise storage
        x = chebpts(nMax1);            % Chebyshev points
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
        % Convert the discrete values to a CHEBFUN.
        p = chebfun(P(:,cc), dom); 
    
    case 2 % QR

        pts = 2*nMax1;              % Expand on Chebyshev grid of twice the size
        [x, w] = chebtech2.chebpts(pts);      % Grab the Clenshaw-Curtis weights
        theta = pi*(pts-1:-1:0)'/(pts-1);
        A = cos(theta*(0:nMax));              % Vandemonde-type matrix
        D = spdiags(sqrt(w(:)), 0, pts, pts); % C-C quad weights
        Dinv = spdiags(1./sqrt(w(:)), 0, pts, pts);
        [Q, R] = qr(D*A, 0);                  % Weighted QR
        P = Dinv*Q;
        if ( normalize )
            PP = P(:,n+1) * diag(sqrt(2/diff(dom)) * sign(P(end,n+1)));
        else
            PP = P(:,n+1) * diag(1./P(end,n+1));
        end
        % Convert discrete values to a CHEBFUN:
        p_chebtech = chebtech2(PP);              % Make a CHEBTECH.
        p_chebtech = prolong(p_chebtech, nMax1); % Force the correct length.
        p_bndfun = bndfun(p_chebtech, dom);      % Make a BNDFUN.
        p = chebfun({p_bndfun});                 % Make a CHEBFUN.
            
    case 3 % LEG2CHEB

        c_leg = [1 ; zeros(n, 1)];              % Legendre coefficients.
        c_cheb = chebtech.leg2cheb(c_leg);      % Chebyshev coefficients.
        if ( normalize )
            c_cheb = c_cheb*sqrt((n+.5));
        end
        p_chebtech = chebtech2({[], c_cheb});   % Make a CHEBTECH.
        p_bndfun = bndfun(p_chebtech, dom);     % Make a BNDFUN.
        p = chebfun({p_bndfun});                % Make a CHEBFUN.

end

% Adjust orientation:
if ( size(n, 1) > 1 )
   p = p.'; 
end

