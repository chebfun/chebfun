function f = compose(f, op, g, data, pref)
%COMPOSE   Composition of CHEBTECH objects.
%   COMPOSE(F, OP) returns a CHEBTECH representing OP(F) where F is also a
%   CHEBTECH object, and OP is a function handle.
%
%   COMPOSE(F, OP, G) returns a CHEBTECH representing OP(F, G) where F and G
%   are CHEBTECH objects, and OP is a function handle.
%
%   COMPOSE(F, G) returns a CHEBTECH representing G(F), where both F and G are
%   also CHEBTECH objects. If the range of F is not in [-1, 1] then an error is
%   thrown.
%
%   COMPOSE(F, OP, G, DATA, PREF) or COMPOSE(F, OP, [], DATA, PREF) uses the
%   constructor data in the structure DATA and the options passed by the
%   preferences structure PREF to build the returned CHEBTECH.  In particular,
%   one can set PREF.REFINEMENTFUNCTION to be a function which takes advantage
%   of F and possibly OP or G being CHEBTECH objects.

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
pref.minSamples = max(pref.minSamples, length(f));
pref.eps = max(pref.eps, f.epslevel);
pref.sampleTest = false;

if ( nfuns == 2 )
    
    if ( size(f, 2) ~= size(g, 2) )
        error('CHEBFUN:CHEBTECH:compose:dim', ...
              'Matrix dimensions must agree.')
    end

    % Grab some data from G:
    pref.minSamples = max(pref.minSamples, length(g));
    pref.eps = max(pref.eps, g.epslevel);
    
elseif ( isa(op, 'chebtech') )
       
    if ( (size(op, 2) > 1) && (size(f, 2) > 1) )
        error('CHEBFUN:CHEBTECH:compose:arrval', ...
              'Cannot compose two array-valued CHEBTECH objects.')
    end

    % If OP is a CHEBTECH, we grab some of its data:
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

% Make CHEBTECH object:
if ( ~isfield(data, 'hscale') || isempty(data.hscale) )
    data.hscale = f.hscale;
end
f = f.make(op, data, pref);

end
