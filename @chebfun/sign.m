function g = sign(f, varargin)
%SIGN   Sign function.
%   G = SIGN(F) returns a piecewise constant chebfun G such that G(x) = 1 in the
%   interval where F(x) > 0, G(x) = -1 in the interval where F(x) < 0 and G(x) =
%   0 in the interval where F(x) = 0. The breakpoints of H are introduced at
%   zeros of F.
%
%   G = SIGN(F, 'reassignimpulses') reassigns the impulses.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Copy g to f:
g = f;

% Deal with the empty case:
if ( isempty(f) )
    return
end

% If f is not real, sign returns f?
if ( ~isreal(f) )
    % [TODO]: Perhaps it should return:
%     g = f./abs(f);
    return
end

% [TODO]: If f is the zero Chebfun, then return zero:
% g.impulses = 0*g.impulses(1,:);
% if ( ~any(g) )
%     % There are nonzero impulses:
%     if ( any(f.impulses) )
%         imps = sign(double(f.impulses(1,:)));
%         idx = [true, logical(imps(2:end-1)), true];
%         g = chebfun(0, f.domain(idx));
%         g.imps = imps(idx);
%     end
%     return
% end

% Set scales and tolerances:
vs = get(f, 'vscale'); vs = max([vs{:}]); % [TODO]
hs = hscale(f);
pref = chebfun.pref();
tol = 100*pref.chebfun.eps*vs;
dom = f.domain;

% Compute the roots of the chebfun:
r = roots(f, 'nozerofun');

% Include the endpoints in the roots list: (since this forms the new domain)
if ( isempty(r) )
    r = dom([1 end]).';
else
    if ( abs(r(1)  - dom(1)) > 1e-14*hs )
        % Leftmost root is far from domain(1):
        r = [ dom(1) ; r  ];
    else
        r(1) = dom(1);
    end
    if ( abs(r(end) - dom(end)) > 1e-14*hs )
        % Rightmost root is far from domain(1):
        r = [ r ; dom(end) ];
    else
        r(end) = dom(end);
    end
end

% Check for double roots: (double roots may be quite far apart!)
ind = find(diff(r) < 1e-7*hs);
cont = 1;
while ( ~isempty(ind) && cont < 3 )
    remove = false(size(r));
    for k = 1:length(ind)
        % Check whether a double root or two single roots very close
        if ( abs(feval(f, mean(r(ind(k):ind(k)+1)))) < tol )
            if ( ismember(r(ind(k)+1), dom) )
                remove(ind) = true;
            else
                remove(ind+1) = true;
            end
        end
    end
    r(remove) = [];
    cont = cont + 1;
    ind = find(diff(r) < 1e-7*hs);
end
% Make sure that the domain of definition is not changed
r([1, end]) = dom([1, end]);

% Non-trivial impulses should not be removed:
if ( nargin == 1 )
    lvals = dom;
    for k = 1:length(dom)-1
        lvals(k) = get(f.funs{k}, 'lval');
    end
    lvals(end) = get(f.funs{end}, 'rval');
    rvals = dom;
    for k = 2:length(dom)
        rvals(k) = get(f.funs{k-1}, 'rval');
    end;
    rvals(1) = get(f.funs{1}, 'lval');
    idxl = abs(f.impulses(1,:) - lvals) > 100*tol;
    idxr = abs(f.impulses(1,:) - rvals) > 100*tol;
    idx = idxl | idxr;
    idx([1 end]) = 1;
    r = sort(unique([r ; dom(idx).']));
    [ignored, idx2] = intersect(r,dom(idx));
end

% Evaluate at an arbitrary point between each root:
c = 0.5912;
nr = length(r);
vals = sign( feval( f , c*r(1:nr-1) + (1-c)*r(2:nr) ) );

% Reshape r as a row vector:
r = r(:).';

% Build new chebfun
ff = cell(1, nr-1);
for k = 1:nr-1
    ff{k} = fun.constructor(vals(k), r(k:k+1));
end
g.funs = ff;
g.domain = r;
imps = f.impulses;
g.impulses = chebfun.jumpVals(ff);

% Reassign impulses
if ( nargin == 1 )
    g.impulses(1,idx2) = sign( imps(1, idx) );
end

end