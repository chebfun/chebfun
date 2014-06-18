function [h, ishappy] = merge(f, g, vscale, hscale, pref)
%MERGE   Attempt to merge to FUN objects into a single object.
%   MERGE(F, G) attempts to merge two FUN objects on neighbouring domains to a
%   single FUN object on the combined domain. A warning is issued if the merge
%   isnot successful.
%
%   [H, ISHAPPY] = MERGE(F,G) prevents the warning being thrown and returns a
%   boolean value for happiness instead.
%
%   ... = MERGE(F, G, VSCALE, HSCALE, PREF) uses uses the values VSCALE and
%   HSCALE and the CHEBFUNPREF object PREF in the construction of the new FUN.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Extract domain info:
domF = f.domain;
domG = g.domain;
dom = [domF, domG];
newDom = dom([1,end]);

% Parse inputs:
if ( (nargin < 3) || isempty(vscale) )
    vscale = 0;
end
if ( (nargin < 4) || isempty(hscale) )    
    hscale = norm(dom(isfinite(dom)), inf);
end
if ( nargin < 5 )
    pref = chebfunpref();
    pref.eps = max(get(f, 'epslevel'), get(g, 'epslevel'));
end
tol = max(pref.eps);

% Check the domains:
if ( abs(domF(2) - domG(1)) > hscale*tol )
    error('CHEBFUN:FUN:merge:domains', ...
        'F and G must be on consecutive domains.');
end

% Assign new domain, vscale, and hscale to a data struct to pass to constructor:
data = struct('domain', newDom, 'vscale', vscale, 'hscale', hscale);
   
% Grab the correct exponents:
if ( issing(f) && issing(g) )
    expsF = get(f, 'exponents');
    expsG = get(g, 'exponents');
    data.exponents = [expsF(1), expsG(2)];
elseif ( issing(f) )
    expsF = get(f, 'exponents');
    data.exponents = [expsF(1), 0];
elseif ( issing(g) )
    expsG = get(g, 'exponents');
    data.exponents = [0, expsG(2)];
end
    
% Attempt to form a merged FUN:
h = fun.constructor(@(x) myFun(x, f, g, dom), data, pref);

ishappy = get(h, 'ishappy');
if ( ~ishappy && (nargout < 2) )
    warning('CHEBFUN:FUN:merge:unsuccessful', ...
        'Attempted merge was unnsuccessful with %d points.', length(h));
end

end

function y = myFun(x, f, g, dom)
%MYFUN(X, F, G, DOM) evaluates F or G, depending in which domain X lies.

% Evaluate F:
idxF = x <= dom(2);
if ( any(idxF) )
    y(idxF,:) = feval(f, x(idxF));
end
% Evaluate G:
idxG = x >= dom(3);
if ( any(idxG) )
    y(idxG,:) = feval(g, x(idxG));
end

end
