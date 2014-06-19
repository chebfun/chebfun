function f = compose(f, op, g, data, pref)
%COMPOSE   Composition of FOURTECH objects.
%   COMPOSE(F, OP) returns a FOURTECH representing OP(F) where F is also a
%   FOURTECH object, and OP is a function handle.
%
%   COMPOSE(F, OP, G) returns a FOURTECH representing OP(F, G) where F and G
%   are FOURTECH objects, and OP is a function handle.
%
%   FOURTECH(F, G) returns a FOURTECH representing G(F), where both F and G are
%   also FOURTECH objects. If the range of F is not in [-1, 1] then an error is
%   thrown.
%
%   COMPOSE(F, OP, G, PREF) or COMPOSE(F, OP, [], PREF) uses the options passed
%   by the preferences structure PREF to build the returned FOURTECH. In
%   particular, one can set PREF.REFINEMENTFUNCTION to be a function which takes
%   advantage of F and possibly OP or G being FOURTECH objects.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs:
if ( nargin < 5 )
    pref = f.techPref();
else
    pref = f.techPref(pref);
end

if ( nargin < 4 )
    data = struct();
end

if ( (nargin < 3) || isempty(g) )
    nfuns = 1;
else
    nfuns = 2;
end

% Set some preferences:
vscale = f.vscale;
pref.minSamples = max(pref.minSamples, length(f));
pref.eps = max(pref.eps, f.epslevel);
pref.sampleTest = false;

if ( nfuns == 2 )
    if ( size(f, 2) ~= size(g, 2) )
        error('CHEBFUN:FOURTECH:compose:dim', ...
              'Matrix dimensions must agree.')
    end

    % Grab some data from G2:
    vscale = max(vscale, g.vscale);
    pref.minSamples = max(pref.minSamples, length(g));
    pref.eps = max(pref.eps, g.epslevel);
    
elseif ( isa(op, 'fourtech') )
    % If OP is a FOURTECH, we grab some of its data:
    if ( (size(op, 2) > 1) && (size(f, 2) > 1) )
        error('CHEBFUN:FOURTECH:compose:arrval', ...
              'Cannot compose two array-valued FOURTECH objects.')
    end

    if ( norm(f.values(:), inf) > 1 )
        error('CHEBFUN:FOURTECH:compose:range', ...
              'The range of f is not contained in the domain of g.')
    end

    vscale = max(vscale, op.vscale);
    pref.minSamples = max(pref.minSamples, length(op));
    pref.eps = max(pref.eps, op.epslevel);
    
end

% Use a naive evaluation procedure if a custom refinement has not been passed.
if ( ischar(pref.refinementFunction) )
    if ( nfuns == 1 )
        op = @(x) feval(op, feval(f, x));
    else
        op = @(x) feval(op, feval(f, x), feval(g, x));
    end
end

% Make FOURTECH object:
if ( ~isfield(data, 'vscale') || isempty(data.vscale) )
    data.vscale = vscale;
end

if ( ~isfield(data, 'hscale') || isempty(data.hscale) )
    data.hscale = f.hscale;
end

f = f.make(op, data, pref);
f = simplify(f);

% Throw a warning:
if ( ~f.epslevel )
    warning('FOURTECH:FOURTECH:compose:convfail', ...
        [ 'Composition with ', func2str(op), ...
          ' failed to converge with ', int2str(length(f)), ' points.' ]);
end

end