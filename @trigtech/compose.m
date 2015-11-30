function f = compose(f, op, g, data, pref)
%COMPOSE   Composition of TRIGTECH objects.
%   COMPOSE(F, OP) returns a TRIGTECH representing OP(F) where F is also a
%   TRIGTECH object, and OP is a function handle.
%
%   COMPOSE(F, OP, G) returns a TRIGTECH representing OP(F, G) where F and G
%   are TRIGTECH objects, and OP is a function handle.
%
%   TRIGTECH(F, G) returns a TRIGTECH representing G(F), where both F and G are
%   also TRIGTECH objects. If the range of F is not in [-1, 1] then an error is
%   thrown.
%
%   COMPOSE(F, OP, G, PREF) or COMPOSE(F, OP, [], PREF) uses the options passed
%   by the preferences structure PREF to build the returned TRIGTECH. In
%   particular, one can set PREF.REFINEMENTFUNCTION to be a function which takes
%   advantage of F and possibly OP or G being TRIGTECH objects.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
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
pref.eps = max(pref.eps, eps);
pref.sampleTest = false;

if ( nfuns == 2 )
    if ( size(f, 2) ~= size(g, 2) )
        error('CHEBFUN:TRIGTECH:compose:dim', ...
              'Matrix dimensions must agree.')
    end

    % Grab some data from G:
    pref.minSamples = max(pref.minSamples, length(g));
    
elseif ( isa(op, 'trigtech') )
    % If OP is a TRIGTECH, we grab some of its data:
    if ( (size(op, 2) > 1) && (size(f, 2) > 1) )
        error('CHEBFUN:TRIGTECH:compose:arrval', ...
              'Cannot compose two array-valued TRIGTECH objects.')
    end

    if ( norm(f.values(:), inf) > 1 + eps )
        error('CHEBFUN:TRIGTECH:compose:range', ...
              'The range of f is not contained in the domain of g.')
    end

    pref.minSamples = max(pref.minSamples, length(op));
    
end

% Use a naive evaluation procedure if a custom refinement has not been passed.
if ( ischar(pref.refinementFunction) )
    if ( nfuns == 1 )
        op = @(x) feval(op, feval(f, x));
    else
        op = @(x) feval(op, feval(f, x), feval(g, x));
    end
end

f = f.make(op, data, pref);

% Throw a warning:
if ( ~f.ishappy )
    warning('TRIGTECH:TRIGTECH:compose:convfail', ...
        [ 'Composition with ', func2str(op), ...
          ' failed to converge with ', int2str(length(f)), ' points.' ]);
end

end
