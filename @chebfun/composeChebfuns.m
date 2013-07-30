function h = composeChebfuns(f, g, pref)
%COMPOSE chebfun composition
%   COMPOSE(F,G) returns the composition of the chebfuns F and G, F(G). The
%   range of G must be in the domain of F.
%
%   Example:
%           f = chebfun(@(x) 1./(1+100*x.^2));
%           g = chebfun(@(x) asin(.99*x)/asin(.99));
%           h = compose(f,g);

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( nargin < 3 )
    pref = chebfun.pref();
else
    pref = chebfun.pref(pref);
end

if ( ~isa(g,'chebfun') || ~isa(f, 'chebfun') )
    warning('CHEBFUN:composeChebfuns:notachebfun', ...
        'Inputs should be chebfun objects!');
    h = compose(f, g, pref);
    return
end

% [TODO]
% % Check to see if f or g is 'x';
% xg = chebfun(@(x) x, g.domain);
% if ( norm(g - xg, 2) < get(g, 'epslevel') )
%     h = f;
%     return
% end
% xf = chebfun(@(x) x, f.domain);
% if ( norm(f - xf, 2) < get(g, 'epslevel') )
%     h = g;
%     return
% end

% Delta functions ?
if ( size(f.impulses, 1) > 1 || size(g.impulses, 1) > 1 )
    warning('CHEBFUN:compose:imps',  ...
        'Composition does not handle delta functions')
end

% % g must be a real-valued function
% if ~isreal(g)
%     %     error('CHEBFUN:compose:complex', 'G must be real valued to construct F(G).')
%     %     warning('CHEBFUN:compose:complex', 'G SHOULD be real valued to construct F(G).');
%     % Experimental feature allows composition when G has a complex range.
%     %   This is only really of any use when F is constructed from a
%     %   polynomial otherwise approximation off the real line is awful.
% end

tol = 100*pref.chebfun.eps;

% Range of g must be in the domain of f.
mm = minandmax(g);
if ( f.domain(1) > mm(1) + tol || f.domain(end) < mm(2) - tol && isreal(g) )
    error('CHEBFUN:compose:domain', ...
        'Range of G, [%g, %g], must be in the domain of F, [%g, %g].', ...
        mm(1), mm(2), f.ends(1), f.ends(2))
end

% If f has breakpoints, find the corresponding x-points in the domain of g:
bkpts = [];
if ( numel(f.domain) > 2 )
    bkf = f.domain(f.domain > mm(1) - tol & f.domain < mm(2) + tol);
    for k = 1:length(bkf)
        bkpts = [bkpts; roots(g - bkf(k))];
    end
end
dom = union(g.domain, bkpts);

% funs = {};
for k = 1:length(dom)-1
%     newfun = composeFuns(f.funs{k}, gg.funs{k},
    gk = restrict(g, dom(k:k+1));
    h = compose(gk, @(g) feval(f, g));
end

% Fix imps values
h.impulses(1,:) = feval(f, feval(g, h.domain));
