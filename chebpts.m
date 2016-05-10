function [x, w, v, t] = chebpts(n, dom, type)
%CHEBPTS    Chebyshev points.
%   CHEBPTS(N) returns N Chebyshev points of the 2nd-kind in [-1,1].
%
%   CHEBPTS(N, D), where D is vector of length 2 and N is a scalar integer,
%   scales the nodes and weights for the interval [D(1),D(2)]. If length(D) > 2
%   and N a vector of length(D)-1, then CHEBPTS(N, D) returns a column vector of
%   the stacked N(k) Chebyshev points on the subintervals D(k:k+1). If length(N)
%   is 1, then D is treated as [D(1),D(end)].
%
%   [X, W] = CHEBPTS(N, D) returns also a row vector of the (scaled) weights for
%   Clenshaw-Curtis quadrature (computed using [1]). (For nodes and weights of
%   Gauss-Chebyshev quadrature, use [X, W] = JACPTS(N, -.5, -.5, D))
%
%   [X, W, V] = CHEBPTS(N, D) returns, in addition to X and W, the barycentric
%   weights V corresponding to the Chebyshev points X. The weights are scaled to
%   have infinity norm 1.
%
%   [X, W, V, T] = CHEBPTS(N) returns also the angles T so that cos(T) = X.
%
%   [X, W, V, T] = CHEBPTS(N, KIND) or CHEBPTS(N, D, KIND) returns Chebyshev
%   points, weights, and angles of the 1st-kind if KIND = 1 and 2nd-kind if KIND
%   = 2 (default).
%
% See also TRIGPTS, LEGPTS, JACPTS, LAGPTS, HERMPTS, LOBPTS, and RADAUPTS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Mathematical reference]:
%   Jarg Waldvogel, "Fast construction of the Fejer and Clenshaw-Curtis
%   quadrature rules", BIT Numerical Mathematics, 46, (2006), pp 195-202.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse inputs:
if ( nargin == 2 )
    if ( length(dom) == 1 )
        type = dom;
        dom = chebfunpref().domain;
    else
        type = 2;
    end
end
if ( nargin == 1 )
    dom = chebfunpref().domain;
    type = 2;
end

% Verify that the number of points requested and the domains match:
if ( length(n) == 1 && length(dom) > 1 )
    % Ignore interior breaks in this instance.
    dom = dom([1,end]); 
elseif ( length(n) ~= length(dom) - 1 )
    error('CHEBFUN:chebpts:mismatchND', 'Vector N does not match domain D.'); 
end

% Create a dummy CHEBTECH of appropriate type to access static methods.
f = feval(['chebtech', num2str(type)]);

if ( length(n) == 1 ) % Single grid.
    
    % Call the static CHEBTECH.CHEBPTS() method:
    x = f.chebpts(n);
    % Scale the domain:
    x = scaleNodes(x, dom);
    
    if ( nargout > 1 )
        % Call the static CHEBTECH.CHEBPTS() method:
        w = f.quadwts(n);
        % Scale the domain:
        w = scaleWeights(w, dom);
    end
    if ( nargout > 2 )
        v = f.barywts(n);
        
        % Rescale the barycentric weights for 1st-kind points so that 
        % norm(v, inf) = 1:
        if ( type == 1 )
            v = v/norm(v, inf);
        end
    end
    if ( nargout > 3 )
        t = f.angles(n);
    end
    
else                  % Piecewise grid.
    
    % Initialise cell arrays for temporary storage:
    x = cell(length(n), 1);
    w = cell(1, length(n));
    v = cell(length(n), 1);
    % Loop over the number of subintervals: (as above)
    for k = 1:numel(n)
        x{k} = f.chebpts(n(k));
        x{k} = scaleNodes(x{k}, dom(k:k+1));
        if ( nargout > 1 )
            w{k} = f.quadwts(n(k));
            w{k} = scaleWeights(w{k}, dom(k:k+1));
        end
        if ( nargout > 2 )
            v{k} = f.barywts(n(k));
        end
    end
    % Convert the cell to an array for output:
    x = cell2mat(x);
    w = cell2mat(w);
    v = cell2mat(v);
    if ( nargout > 3 )
        a = dom(1);
        b = dom(end);
        tmp = (x - a)/(b - a) - (b - x)/(b - a);
        t = acos(tmp);
    end
    
end

end

function y = scaleNodes(x, dom)
%SCALENODES   Scale the Chebyshev nodes X from [-1,1] to DOM.
% TODO: Deal with unbounded domains
if ( dom(1) == -1 && dom(2) == 1 )
    % Nodes are already on [-1, 1];
    y = x;
    return
end
% Scale the nodes:
y = dom(2)*(x + 1)/2 + dom(1)*(1 - x)/2;
end

function w = scaleWeights(w, dom)
%SCALEWEIGHTS   Scale the Chebyshev weights W from [-1,1] to DOM.
% TODO: Deal with unbounded domains
if ( dom(1) == -1 && dom(2) == 1 )
    % Nodes are already on [-1, 1];
    return
end
% Scale the weights:
w = (diff(dom)/2)*w;
end
