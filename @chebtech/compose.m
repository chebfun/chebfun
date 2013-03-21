function f = compose(f, op, g, pref)
%COMPOSE   Composition of CHEBTECH objects.
%   COMPOSE(F, OP) returns a CHEBTECH representing OP(F) where F is also a
%   CHEBTECH2 object, and OP is a function handle.
%
%   COMPOSE(F, OP, G) returns OP(F, G) where F and G are CHEBTECH objects,
%   and OP is a function handle.
%
%   COMPOSE(F, G) returns a CHEBTECH representing G(F), where both F and G are
%   also CHEBTECH objects. If the range of F is not in [-1, 1] then an error is
%   thrown.
%
%   COMPOSE(F, OP, PREF), COMPOSE(F, OP, G, PREF), or COMPOSE(F, G, PREF) uses
%   the options passed by the preferences structure PREF. In particular, one can
%   pass a PREF.CHEBTECH.REFINMENTFUNCTION which takes advantage of the fact
%   that F (and possibly OP or G) are CHEBTECH objects.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

nfuns = 2;
% Parse inputs:
if ( (nargin > 2) && isstruct(g) )
    pref = g;
    g = [];
    nfuns = 1;
elseif ( nargin < 4 )
    pref = f.pref();
end
if ( (nargin < 3) || isempty(g) )
    nfuns = 1;
    g = [];
end

% Set some preferences:
vscale = f.vscale;
pref.chebtech.minSamples = max(pref.chebtech.minSamples, length(f));
pref.chebtech.eps = max(pref.chebtech.eps, f.epslevel);
pref.chebtech.sampletest = false;
if ( nfuns == 2 )
    if ( size(f, 2) ~= size(g, 2) )
        error('CHEBFUN:CHEBTECH:compose:dim', ...
            'Matrix dimensions must agree.')
    end
    % Grab some data from G2:
    vscale = max(vscale, g.vscale);
    pref.chebtech.minSamples = max(pref.chebtech.minSamples, length(g));
    pref.chebtech.eps = max(pref.chebtech.eps, g.epslevel);
elseif ( isa(op, 'chebtech') )
    % If OP is a CHEBTECH, we grab some of its data:
    if ( (size(op, 2) > 1) && (size(f, 2) > 1) )
        error('CHEBFUN:CHEBTECH:compose:multival', ...
            'Cannot compose two multivalued CHEBTECH objects.')
    end
    if ( norm(f.values(:), inf) > 1 )
        error('CHEBFUN:CHEBTECH:compose:range', ...
            [ 'The range of f (approx [' num2str(min(f.values)), ', ', ...
            num2str(max(f.values)), ']) is not in the domain of G ([-1,1])' ])
    end
    vscale = max(vscale, op.vscale);
    pref.chebtech.minSamples = max(pref.chebtech.minSamples, length(op));
    pref.chebtech.eps = max(pref.chebtech.eps, op.epslevel);
%     op = @(x) feval(op, x); % This isn't needed, as we use FEVAL in refFunc.
% [TODO]: remove line above?
end

% Use a naive evaluation procedure if a custon refinement has not been passed.
if ( ischar(pref.chebtech.refinementFunction) )
    if ( nfuns == 1 )
        op = @(x) feval(op, feval(f, x));
    else
        op = @(x) feval(op, feval(f, x), feval(g, x));
    end   
end

% Make CHEBTECH object:
f = f.make(op, vscale, f.hscale, pref);

% Throw a warning:
if ( ~f.epslevel )
    warning('CHEBFUN:CHEBTECH:compose:convfail', ...
        [ 'Composition with ', func2str(op), ...
        ' failed to converge with ', int2str(length(f)), ' points.' ]);
end

end
