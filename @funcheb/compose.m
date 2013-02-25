function f = compose(f, op, g, pref)
%COMPOSE  Compostition of FUNCHEB objects.
%   COMPOSE(F, OP) returns a FUNCHEB representing OP(F) where F is also a
%   FUNCHEB2 object, and OP is a function handle.
%
%   COMPOSE(F, OP, G) returns OP(F, G) where F and G are FUNCHEB objects,
%   and OP is a function handle.
%
%   COMPOSE(F, G) returns a FUNCHEB representing G(F), where both F and G are
%   also FUNCHEB objects. If the range of F is not in [-1, 1] then an error is
%   thrown.
%
%   COMPOSE(F, OP, PREF), COMPOSE(F, OP, G, PREF), or COMPOSE(F, G, PREF) uses
%   the options passed by the prefences structure PREF. In particular, one can
%   pass a PREF.(class(F)).refinmentFunction which takes advantage of the fact
%   that F (and possibly OP or G) are FUNCHEB objects.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

nfuns = 2;
% Parse inputs:
if ( nargin > 2 && isstruct(g) )
    pref = g;
    g = [];
    nfuns = 1;
elseif ( nargin < 4 )
    pref = f.pref();
end
if ( nargin < 3 || isempty(g) )
    nfuns = 1;
    g = [];
end

% Set some preferences:
vscale = f.vscale;
pref.(class(f)).minSamples = max(pref.(class(f)).minSamples, length(f));
pref.(class(f)).eps = max(pref.(class(f)).eps, f.epslevel);
pref.(class(f)).sampletest = false;
if ( nfuns == 2 )
    if ( size(f, 2) ~= size(g, 2) )
        error('CHEBFUN:FUNCHEB:compose:dim', ...
            'Matrix dimensions must agree.')
    end
    % Grab some data from G2:
    vscale = max(vscale, g.vscale);
    pref.(class(f)).minSamples = max(pref.(class(f)).minSamples, length(g));
    pref.(class(f)).eps = max(pref.(class(f)).eps, g.epslevel);
elseif ( isa(op, 'funcheb') )
    % If OP is a FUNCHEB, we grab some of its data:
    if ( size(op, 2) > 1 && size(f, 2) > 1)
        error('CHEBFUN:FUNCHEB:compose:multival', ...
            'Cannot compose two multivalued FUNCHEB2 objects.')
    end
    if ( norm(f.values(:), inf) > 1 )
        error('CHEBFUN:FUNCHEB:compose:range', ...
            ['The range of f (approx [' num2str(min(f.values)), ', ', ...
            num2str(max(f.values)), ']) is not in the domain of G ([-1,1])'])
    end
    vscale = max(vscale, op.vscale);
    pref.(class(f)).minSamples = max(pref.(class(f)).minSamples, length(op));
    pref.(class(f)).eps = max(pref.(class(f)).eps, op.epslevel);
%     op = @(x) feval(op, x); % This isn't needed, as we use FEVAL in refFunc.
end

% Use a naive evaluation procedure is a custon refinement has not been passed.
if ( ischar(pref.(class(f)).refinementFunction) )
    if ( nfuns == 1 )
        op = @(x) feval(op, feval(f, x));
    else
        op = @(x) feval(op, feval(f, x), feval(g, x));
    end   
end

% Call the funcheb2 constructor:
f = f.make(op, vscale, f.hscale, pref);

% Throw a warning:
if ( ~f.epslevel )
    warning('CHEBFUN:FUNCHEB:compose:convfail',...
        ['Composition with ', func2str(op), ...
        ' failed to converge with ', int2str(length(f)), ' points.']);
end

end
