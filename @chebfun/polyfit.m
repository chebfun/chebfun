function f = polyfit(y, n, varargin)
%POLYFIT   Fit polynomial to a CHEBFUN.
%   F = POLYFIT(Y, N) returns a CHEBFUN F corresponding to the polynomial of
%   degree N that fits the CHEBFUN Y in the least-squares sense.
%
%   If y is a global polynomial of degree n then this code has an O(n (log n)^2)
%   complexity. If y is piecewise polynomial then it has an O(n^2) complexity.
%
%   F = POLYFIT(CHEBFUN, X, Y, N, D) returns a CHEBFUN F on the domain D which
%   corresponds to the polynomial of degree N that fits the data (X, Y) in the
%   least-squares sense. If D is not given, it is assumed to be [X(1), X(end)].
%
%   Note CHEBFUN/POLYFIT does not not support more than one output argument in
%   the way that MATLAB/POLYFIT does.
%
% See also INTERP1.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% TODO: Come up with a better way of accessing POLYFIT(X, Y, N, D).

if ( nargout > 1 )
    error('CHEBFUN:polyfit:nargout', ...
        'Chebfun/polyfit only supports one output');
end

if ( nargin > 2 )
    % POLYFIT(CHEBFUN, X, Y, N, D)
    f = discretePolyfit(n, varargin{:});
    return
else
    % POLYFIT(Y, D)
    f = continuousPolyfit(y, n);
    return
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = discretePolyfit(x, y, n, dom)
%DISCRETEPOLYFIT
% [TODO]: Document.
if ( nargin < 4 )
    % Use x data to determine domain:
    dom = [x(1), x(end)];
end

% % Call built in POLYFIT():
% p = polyfit(x, y, n);
% % Create a CHEBFUN:
% f = chebfun(@(x) polyval(p, x), dom, length(p));

% Use BARY():
w = baryWeights(x);                    % Barycentric weights for x.
xCheb = chebpts(n);                    % Chebyshev grid.
yCheb = chebtech.bary(xCheb, y, x, w); % Evaluate interpolant at xCheb.
f = chebfun(yCheb, dom);               % Create a CHEBFUN of the result.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = continuousPolyfit(y, n)
%CONTINUOUSPOLYFIT
% [TODO]: Document.

% The code below is a fast version of the code:
% Enorm = legpoly(0:n,[a,b],'norm');  % Legendre-Vandermonde matrix
% f = Enorm*(Enorm'*y);               % least squares chebfun

if ( any(isinf(y.domain)) )
    error('CHEBFUN:polyfit:unbounded', 'Unbounded domains are not supported.');
end

% Copy y to f:
f = y;

if ( n > length(y) && numel(y.funs) == 1 )
    % Nothing to do here!
    return
end

if ( min(size(y)) > 1 )
    error('CHEBFUN:polyfit:arrayvalues', ...
        'Array-valued CHEBFUN objects are not supported by POLYFIT().');
end

if ( numel(y.funs) == 1 )
    f.funs{1} = polyfit(y.funs{1}, n);
    return
end

a = y.domain(1);
b = y.domain(end);

isSimple = all(cellfun(@(f) isa(f.onefun, 'chebtech'), y.funs));
if ( isSimple )
    
    scl = 2./(2*(0:n)'+1);  % Orthonormal scaling.
    cleg = zeros(n+1, 1);   % Initialise storage.
    
    % For each subinterval calculate int P_k f(x)dx and add them up:
    for j = 1:numel(y.funs)
        
        yfun = y.funs{j};   % jth fun.
        xdom = yfun.domain; % Subinterval j:
        % We must map LEGPTS() from [-1, 1] to [a, b] to [xdom].
        zdom = 2*(xdom - a)/(b - a) - 1;
        % Gauss-Legendre nodes and weights on scaled subinterval (zdom):
        [z, w] = legpts(max(length(yfun), n+1), zdom);
        % Mapping from zdom to xdom:
        x = (z+1)*(b-a)/2 + a;
        % Evaluate on y on xdom:
        vals = feval(yfun, x);
        
        % Compute the first two innerproducts by hand:
        cleg(1,:) = cleg(1,:) + w*vals;
        cleg(2,:) = cleg(2,:) + w*(z.*vals);
        % Evaluate Legendre-Vandermonde matrix by recureence relation:
        Pm2 = 1; Pm1 = z;
        for kk = 1:n-1
            P = (2-1/(kk+1))*Pm1.*z - (1-1/(kk+1))*Pm2;
            Pm2 = Pm1; Pm1 = P;
            % Add contribution from subinterval to (k+1)st coefficient:
            cleg(kk+2) = cleg(kk+2) + w*(P.*vals);
        end
        
    end
    
    % Convert the computed Legendre expansion to a Chebysev series:
    cleg = flipud(cleg ./ scl);                  % Scale appropriately.
    c = chebtech.leg2cheb(cleg);                 % Compute Chebyshev coeffs.
    f_chebtech = y.funs{1}.onefun.make({[], c}); % Make a CHEBTECH.
    f_fun = fun.constructor(f_chebtech, [a, b]); % Make a BNDFUN.
    f = chebfun({f_fun});                        % Make a CHEBFUN.
%     f = chebfun(c, 'coeffs', [a, b]);            % Make a CHEBFUN.
    
else
    
    % Legendre-Vandermonde matrix:
    Enorm = legpoly(0:n, [a, b], 'norm');
    % Take project and expand in the Legendre basis:
    f = Enorm*(Enorm'*y);
    
end

end
