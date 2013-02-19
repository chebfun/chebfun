function f = compose(f, op, g, pref)
%COMPOSE  Compostition of FUNCHEB1 objects.
%   COMPOSE(F, OP) returns a FUNCHEB1 representing OP(G) where G is also a
%   FUNCHEB1 object, and OP is a function handle.
%
%   COMPOSE(G1, OP, G2) returns OP(G1, G2) where F and G are FUNCHEB1 objects,
%   and OP is a function handle.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

nfuns = 2;

% Parse inputs
if ( nargin > 2 && isstruct(g) )
    pref = g;
    g = [];
    nfuns = 1;
elseif ( nargin < 4 )
    pref = funcheb1.pref;
end
if ( nargin < 3 || isempty(g) )
    nfuns = 1;
    g = [];
end

% Set some preferences:
vscale = f.vscale;
pref.funcheb1.minSamples = max(pref.funcheb1.minSamples, length(f));
pref.funcheb1.eps = max(pref.funcheb1.eps, f.epslevel);
pref.funcheb1.sampletest = false;
if ( nfuns == 2 )
    if ( size(f, 2) ~= size(g, 2) )
        error('CHEBFUN:FUNCHEB1:compose:dim', ...
            'Matrix dimensions must agree.')
    end
    % Grab some data from G2:
    vscale = max(vscale, g.vscale);
    pref.funcheb1.minSamples = max(pref.funcheb1.minSamples, length(g));
    pref.funcheb1.eps = max(pref.funcheb1.eps, g.epslevel);
elseif ( isa(op, 'funcheb1') )
    % If OP is a FUNCHEB1, we grab some of its data:
    if ( size(op, 2) > 1 && size(f, 2) > 1)
        error('CHEBFUN:FUNCHEB1:compose:multival', ...
            'Cannot compose two multivalued FUNCHEB1 objects.')
    end
    vscale = max(vscale, op.vscale);
    pref.funcheb1.minSamples = max(pref.funcheb1.minSamples, length(op));
    pref.funcheb1.eps = max(pref.funcheb1.eps, op.epslevel);
%     op = @(x) feval(op, x); % This isn't needed, as we use FEVAL in refFunc.
end

% Choose t% Choose the appropriate refinement function:
if ( nfuns == 1 )
    refFunc = composeRefFuncSelect(pref, f);
else
    refFunc = composeRefFuncSelect(pref, f, g);
end
pref.funcheb1.refinementFunction = refFunc;

% Call the funcheb1 constructor:
f = funcheb1(op, vscale, pref);

% Throw a warning:
if ( ~f.epslevel )
    warning('CHEBFUN:FUNCHEB1:compose:convfail',...
        ['Composition with ', func2str(op), ...
        ' failed to converge with ', int2str(length(f)), ' points.']);
end

end

function refFunc = composeRefFuncSelect(pref, f, g)
%COMPOSEREFFUNCSELECT  Choose a refinement function
%   COMPOSEREFFUNCSELECT(PREF, F) chooses a refinement function for 
%   constructing OP(F).
%
%   COMPOSEREFFUNCSELECT(PREF, F, G) chooses a refinement function for 
%   constructing OP(F, G).

% Note: Resampling is enforced for Chebyshev points of 1st kind!

if ( nargin == 2 )
    % OP(G1)
    refFunc = @(op, values, pref) composeRefDouble1(op, values, pref, f);
else
    % OP(G1, G2)
    refFunc = @(op, values, pref) composeRefDouble2(op, values, pref, f, g);
end

end

function [values, giveUp] = composeRefDouble1(op, values, pref, f)
%COMPOSEREFDOUBLE1 Refinement function for composing a FUNCHEB1 with resampling.
    
    if ( isempty(values) )
        % Choose initial n based upon minSamples.
        n = 2^ceil(log2(pref.funcheb1.minSamples - 1)) + 1;
    else
        % (Approximately) powers of sqrt(2):
        pow = log2(size(values, 1) - 1);
        if ( pow == floor(pow) && pow > 5 )
            n = round(2^(floor(pow) + .5)) + 1;
            n = n - mod(n, 2) + 1;
        else
            n = 2^(floor(pow)+1) + 1;
        end
    end
    
    % n is too large!
    if ( n > pref.funcheb1.maxSamples )
        giveUp = true;
        return
    else
        giveUp = false;
    end
    
    % Update f values:
    f = prolong(f, n);
    v1 = f.values;
    values = feval(op, v1);

end

function [values, giveUp] = composeRefDouble2(op, values, pref, f, g)
%COMPOSEREFDOUBLE2  Refinement function for composing FUNCHEB1s with resampling.
    
    if ( isempty(values) )
        % Choose initial n based upon minSamples.
        n = 2^ceil(log2(pref.funcheb1.minSamples - 1)) + 1;
    else
        % (Approximately) powers of sqrt(2):
        pow = log2(size(values, 1) - 1);
        if ( pow == floor(pow) && pow > 5 )
            n = round(2^(floor(pow) + .5)) + 1;
            n = n - mod(n, 2) + 1;
        else
            n = 2^(floor(pow)+1) + 1;
        end
    end
    
    % n is too large!
    if ( n > pref.funcheb1.maxSamples )
        giveUp = true;
        return
    else
        giveUp = false;
    end
        
    % Update f and g values:
    f = prolong(f, n);
    v1 = f.values;
    g = prolong(g, n);
    v2 = g.values;
    values = feval(op, v1, v2);

end
