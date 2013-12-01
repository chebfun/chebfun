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
% TODO: This requires testing.

if ( nargout > 1 )
    error('CHEBFUN:polyfit:nargout', ...
        'Chebfun/polyfit only supports one output argument.');
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
% Fit discrete data {X, Y} using a degree N polynomial on the domain DOM.

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
xCheb = chebpts(n, dom);               % Chebyshev grid.
yCheb = chebtech.bary(xCheb, y, x, w); % Evaluate interpolant at xCheb.
f = chebfun(yCheb, dom);               % Create a CHEBFUN of the result.

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = continuousPolyfit(f, n)
%CONTINUOUSPOLYFIT
% Compute best L2 polynomial approximation to a given CHEBFUN F.

if ( any(isinf(f.domain)) )
    error('CHEBFUN:polyfit:unbounded', 'Unbounded domains are not supported.');
end

if ( n > length(f) && numel(f.funs) == 1 && isa(f.funs.onefun, 'chebtech') )
    % Nothing to do here!
    p = f;
    return
end

% Compute first n Legendre coeffs:
cleg = legpoly(f, n).';                      

% Convert to Chebyshev coeffs:
c = zeros(size(cleg));
for k = 1:size(c, 2)
    c(:,k) = chebtech.leg2cheb(cleg(:,k));   
end

% Make a CHEBTECH:
if ( isa(f.funs{1}.onefun, 'chebtech') )
    % Make the CHEBTECH be of the same type as that in F:
    p_chebtech = f.funs{1}.onefun.make({[], c}); 
else
    p_chebtech = chebtech.constructor({[], c});
end

% Make a BNDFUN:
p_fun = fun.constructor(p_chebtech, f.domain([1, end]));    

% Make a CHEBFUN:
p = chebfun({p_fun});                        

end
