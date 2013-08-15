function g = sign(f, varargin)
%SIGN   Sign function of a CHEBFUN.
%   G = SIGN(F) returns a piecewise constant CHEBFUN G such that G(x) = 1 in the
%   interval where F(x) > 0, G(x) = -1 in the interval where F(x) < 0 and G(x) =
%   0 in the interval where F(x) = 0. Breakpoints in G are introduced at zeros
%   of F.
%
%   G = SIGN(F, 'reassignImpulses') reassigns the impulses.
%
% See also ABS, ROOTS.

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
vs = vscale(f);
hs = hscale(f);
el = epslevel(f);
tol = 100*el*vs;
dom = f.domain;
imps = f.impulses;
numFuns = numel(f.funs);
numCols = size(f.funs{1}, 2);

% Compute the roots of the CHEBFUN:
r = roots(f, 'nozerofun');

% Since each column of an array-valued CHEBFUN must have the same breakpoints,
% we simply take unique(r(:)) and remove any remaining NaNs.
r = unique(r(:));
r(isnan(r)) = [];

% % Also discard any roots that are closer than the accuracy of the CHEBFUN:
r([false ; diff(r) < el*hs*vs]) = [];

% Include the endpoints in the roots list: (since this forms the new domain)
% if ( isempty(r) )
%     g = signNoRoots(f);
%     return
% end

% Add the domain boundaries to the roots vector:
if ( ~isempty(r) )
    if ( abs(r(1) - dom(1)) > 1e-14*hs )
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
else
    r = dom([1, end]).';
end

% % [TODO]: Add this back in. prevents introducing a breakpoint in, say, x.^2.
% % Check for double roots (double roots may be quite far apart!):
% ind = find(diff(r) < 1e-7*hs)
% count = 1;
% while ( ~isempty(ind) && count < 3 )
%     remove = false(size(r));
%     for k = 1:length(ind)
%         % Check whether a double root or two single roots very close:
%         x = mean(r(ind(k):ind(k)+1));
%         if ( abs(feval(f, x)) < tol )
%             if ( ismember(r(ind(k)+1), dom) )
%                 remove(ind) = true;
%             else
%                 remove(ind+1) = true;
%             end
%         end
%     end
%     r(remove) = [];
%     count = count + 1;
%     ind = find(diff(r) < 1e-7*hs);
% end

% Make sure that the domain of definition is not changed
r([1 ; end]) = dom([1, end]);

% Non-trivial impulses should not be removed (unless flag is present):
if ( nargin == 1 )
    lvals = zeros(length(domain), numCols);
    for k = 1:numFuns
        lvals(k,:) = get(f.funs{k}, 'lval');
    end
    lvals(numFuns+1,:) = get(f.funs{end}, 'rval');
    rvals = zeros(length(domain), size(f.funs{1}, 2));
    for k = 2:(numFuns+1)
        rvals(k,:) = get(f.funs{k-1}, 'rval');
    end
    rvals(1,:) = get(f.funs{1}, 'lval');
    idxl = abs(imps(:,:,1) - lvals) > 100*tol;
    idxr = abs(imps(:,:,1) - rvals) > 100*tol;
    idx = any(idxl | idxr, 2);
    idx([1 end]) = 1;
    r = sort(unique([r ; dom(idx).']));
    [ignored, idx2] = intersect(r, dom(idx));
end

% Evaluate at an arbitrary point between each root:
c = 0.5912;
nr = length(r);
x = c*r(1:nr-1) + (1-c)*r(2:nr);
vals = sign( feval(f, x) );

% Reshape r as a row vector:
r = r(:).';

% Build new CHEBFUN:
ff = cell(1, nr-1);
for k = 1:nr-1
    ff{k} = fun.constructor(vals(k,:), r(k:k+1));
end
g.funs = ff;
g.domain = r;
g.impulses = chebfun.jumpVals(ff);

% Reassign impulses
if ( nargin == 1 )
    g.impulses(idx2,:) = sign(imps(idx,:,1));
end

end
