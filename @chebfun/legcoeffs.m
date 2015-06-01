function out = legcoeffs(f, varargin)
%LEGCOEFFS    Legendre polynomial coefficients of a CHEBFUN.
%   A = LEGCOEFFS(F, N) returns the first N coefficients in the Legendre series
%   expansion of the CHEBFUN F, so that such that F approximately equals 
%   A(1) P_0(x) + ... + A(N) P_(N-1)(x), where P_N(x) denotes the N-th Legendre
%   polynomial. A is a column vector.
%
%   If F is smooth (i.e., numel(f.funs) == 1), then A = LEGCOEFFS(F) will assume
%   that N = length(F) - 1;
%
%   LEGCOEFFS does not support quasimatrices.
%
% See also CHEBCOEFFS, JACCOEFFS, FOURCOEFFS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( numel(f) > 1 )
    error('CHEBFUN:CHEBFUN:legcoeffs:quasi', ...
        'LEGCOEFFS does not support quasimatrices.');
end

% Call FUN/LEGCOEFFS():
if ( numel( f.funs ) == 1 )
    out = legcoeffs(f.funs{1}, varargin{:});
else
    if ( nargin < 2 )
        error('CHEBFUN:CHEBFUN:legcoeffs:n', ...
            'Input is piecewise, so LEGCOEFFS() expects a second input argument.')
    end
    out = legcoeffsPiecewise(f, varargin{:});
end

end

function out = legcoeffsPiecewise(f, n)
%LEGCOEFFSPIECEWISE  Compute Legendre coefficients of a piecewise smooth CHEBFUN
% If F is 'simple' (i.e., is a piecewise smooth Chebyshev representation), then
% each of the required inner-products are computed so that c_k = int P_k f(x)dx.
%
% If F is non-trivial (e.g., contains exponents or a non-linear map), then the
% coefficients are computed by creating a LEGCOEFFS quasimatrix and computing
% the inner-products using CHEBFUN/MTIMES(). This can be significatnly slower
% than the above.

% Domain:
a = f.domain(1);
b = f.domain(end);

% Orthonormal scaling:
scl = 2./(2*(0:n-1)'+1);               

isSimple = all(cellfun(@(f) isa(f.onefun, 'chebtech'), f.funs));

if ( isSimple )
    % Compute inner-products manually:
    
    out = zeros(n, numColumns(f)); % Initialise storage.
    
    % For each subinterval calculate int P_k f(x)dx and add them up:
    for j = 1:numel(f.funs)
        
        yfun = f.funs{j};   % jth fun.
        xdom = yfun.domain; % Subinterval j:
        % We must map LEGPTS() from [-1, 1] to [a, b] to [xdom].
        zdom = 2*(xdom - a)/(b - a) - 1;
        % Gauss-Legendre nodes and weights on scaled subinterval (zdom):
        [z, w] = legpts(max(length(yfun), n), zdom);
        % Mapping from zdom to xdom:
        x = (z+1)*(b-a)/2 + a;
        % Evaluate on y on xdom:
        vals = feval(yfun, x);
        
        % Compute the first two inner-products by hand:
        out(1,:) = out(1,:) + sum(diag(w)*vals);
        out(2,:) = out(2,:) + sum(diag(w.*z.')*vals);
        % Evaluate Legendre-Vandermonde matrix by recurrence relation:
        Pm2 = 1; Pm1 = z;
        for kk = 1:n-2
            P = (2-1/(kk+1))*Pm1.*z - (1-1/(kk+1))*Pm2;
            Pm2 = Pm1; Pm1 = P;
            % Add contribution from subinterval to (k+1)st coefficient:
            if ( size(vals, 2) > 1 )
                out(kk+2,:) = out(kk+2,:) + ...
                    sum(bsxfun(@(u,v) u*v, w.*P.', vals.').');
            else
                out(kk+2,:) = out(kk+2,:) + (w.*P.')*vals;
            end
        end
        
    end
    
else
    % Compute using LEGCOEFFS(): (This will be much slower!)
    
    % Legendre-Vandermonde matrix:
    Enorm = legpoly(0:n-1, [a, b]);
    % Compute the projection (with correct scaling):
    out = (Enorm.'*f);
    
end

% Scale appropriately:
out = (bsxfun(@rdivide, out, scl));    
    
end
                
