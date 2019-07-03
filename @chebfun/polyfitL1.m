function p = polyfitL1(f, n)
%POLYFITL1 Best polynomial approximation in the L1-norm for real functions.
%
% P = POLYFITL1(F, N) computes the best polynomial approximation to the real
% continuous function F in the L1-sense, using Watson's algorithm. F and P
% are both CHEBFUN objects.
%
% Examples:
%   x = chebfun('x'); f = abs(x);
%   p = polyfitL1(f, 20); plot(f-p)
%
% References:
%
%   [1] Watson, G. A. "An algorithm for linear L1 approximation of
%       continuous functions." IMA Journal of Numerical Analysis,
%       1.2 (1981): 157-167.
%
%   [2] Yuji Nakatsukasa and Alex Townsend. "Error localization of best
%       L1 polynomial approximants." arXiv preprint arXiv:1902.02664 (2019).
%
% See also MINIMAX, CHEBFUN.POLYFIT.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER'S NOTE: Watson's algorithm is a damped Newton iteration that
% (hopefully) converges to the best L1 polynomial approximation to f.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If f is a polynomial and n > deg(f), then p = f:
if ( length( f ) <= n+1 &&  numel( domain(f) ) == 2 )
    p = f;
    return
end

% Fix a tolerance to converge to. Don't be too ambiguous:
tol = 1e-10;

% Grab domain of function:
[a, b] = domain( f );

% Compute the polynomial interpolant pn of f at n+1 Chebyshev points. If
% f-pn has roots only at the Chebyshev points, then Theorem 14.4 and 14.5
% of Powell's Approximation Theory book shows that pn is the best
% polynomial approximation to f in the L1-norm.
x = chebpts(n+3, [a, b]); x = x(2:end-1);
%[TODO] Compute the Chebyshev interpolant using a DCT for speed:
p_interp = chebfun.interp1(x, feval(f, x), [a b]);
r = roots( f - p_interp );

if ( numel( r ) == numel( x ) )
    
    % Very easy case, we already have the best L1 polynomial approximant:
    p = p_interp;
    
else
    
    r = roots(f-p_interp); 
    sums = sign(feval(f-p_interp, a)) * intpsign(n, r, [a,b]);
    sumshistory = norm( sums, inf);
    
    itwatson = 100;
    for ii = 1:itwatson % Watson update
        r = roots(f-p_interp);
        g = sign(feval(f-p_interp, a)) * intpsign(n, r, [a,b]);
        
        A = clenshaw_vec_chebU( 2*(r-a)/(b-a) - 1, eye(n+1)); 
        dfp = diff( f-p_interp );
        D = diag( 2./abs(feval(dfp, r)) );
        H = A'*D*A;
        dc = H \ (g');
        
        % Newton update poly:
        dp = chebfun(ultra2ultra( dc, 1, 0), [a, b], 'coeffs');
        
        % Damped Newton: 
        gam = 1;
        pnow = p_interp;
        while ( gam > 1e-5 )
            p_interp = pnow + gam*dp;
            r = roots( f-p_interp );
            sums = intpsign(n, r, [a,b]);
            if ( norm(sums,inf) < sumshistory && length(r) >= n+1 )
                break
            end
            gam = gam/2;
        end
        sumshistory = norm(sums, inf);
        
        if ( sumshistory < tol / (b-a) )
            break
        end
    end
    if ( ii == itwatson) 
        warning('CHEBFUN:POLYFIT:MAXITER', ... 
                'The maximum number of iterations was reach. Answer may not be accurate.');
    end
    p = p_interp;
end
end

function sums = intpsign(n, r, dom)
%INTPSIGN(N, F, P)

% Grab domain of function:
a = dom(1); b = dom(2); 

% First interval:
[x, w] = legpts( n, [a r(1)] );
xnew = 2*(x-a)/(b-a) - 1;
sums = w*clenshaw_vec_chebU(xnew, eye(n+1));

for jj = 1:length(r)-1
    [x, w] = legpts( n, [r(jj) r(jj+1)] );
    xnew = 2*(x-a)/(b-a) - 1;
    sums = sums + ((-1)^jj)*(w*clenshaw_vec_chebU(xnew, eye(n+1)));
end
if ( length(r) > 1 )
    jj = jj + 1;
    [x, w] = legpts( n, [r(jj) b] );
    xnew = 2*(x-a)/(b-a) - 1;
    sums = sums + ((-1)^jj)*(w*clenshaw_vec_chebU(xnew, eye(n+1)));
end
end

function y = clenshaw_vec_chebU(x, c)
% Clenshaw scheme for evaluating ChebU polynomials: 
x = repmat(x(:), 1, size(c, 2));
bk1 = zeros(size(x, 1), size(c, 2)); 
bk2 = bk1;
e = ones(size(x, 1), 1);
x = 2*x;
n = size(c, 1)-1;
for k = (n+1):-2:3
    bk2 = e*c(k,:) + x.*bk1 - bk2;
    bk1 = e*c(k-1,:) + x.*bk2 - bk1;
end
if ( mod(n, 2) )
    [bk1, bk2] = deal(e*c(2,:) + x.*bk1 - bk2, bk1);
end
y = e*c(1,:) + x.*bk1 - bk2;
end