function f = compose(f, op, g, pref)
%COMPOSE  Compostition of FUNCHEB2 objects.
%   COMPOSE(F, OP) returns a FUNCHEB2 representing OP(F) where F is also a
%   FUNCHEB2 object, and OP is a function handle.
%
%   COMPOSE(F, OP, G) returns OP(F, G) where F and G are FUNCHEB2 objects,
%   and OP is a function handle.
%
%   COMPOSE(F, G) returns a FUNCHEB2 representing G(F), where both F and G are
%   also FUNCHEB2 objects. If the range of F is not in [-1, 1] then an error is
%   thrown.
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
    pref = funcheb2.pref;
end
if ( nargin < 3 || isempty(g) )
    nfuns = 1;
    g = [];
end

% Set some preferences:
vscale = f.vscale;
pref.funcheb2.minSamples = max(pref.funcheb2.minSamples, length(f));
pref.funcheb2.eps = max(pref.funcheb2.eps, f.epslevel);
pref.funcheb2.sampletest = false;
if ( nfuns == 2 )
    if ( size(f, 2) ~= size(g, 2) )
        error('CHEBFUN:FUNCHEB2:compose:dim', ...
            'Matrix dimensions must agree.')
    end
    % Grab some data from G2:
    vscale = max(vscale, g.vscale);
    pref.funcheb2.minSamples = max(pref.funcheb2.minSamples, length(g));
    pref.funcheb2.eps = max(pref.funcheb2.eps, g.epslevel);
elseif ( isa(op, 'funcheb2') )
    % If OP is a FUNCHEB2, we grab some of its data:
    if ( size(op, 2) > 1 && size(f, 2) > 1)
        error('CHEBFUN:FUNCHEB2:compose:multival', ...
            'Cannot compose two multivalued FUNCHEB2 objects.')
    end
    if ( norm(f.values, inf) > 1 )
        error('CHEBFUN:FUNCHEB2:compose:range', ...
            ['The range of f (approx [' num2str(min(f.values)), ', ', ...
            num2str(max(f.values)), ']) is not in the domain of G ([-1,1])'])
    end
    vscale = max(vscale, op.vscale);
    pref.funcheb2.minSamples = max(pref.funcheb2.minSamples, length(op));
    pref.funcheb2.eps = max(pref.funcheb2.eps, op.epslevel);
%     op = @(x) feval(op, x); % This isn't needed, as we use FEVAL in refFunc.
end

% Choose t% Choose the appropriate refinement function:
if ( nfuns == 1 )
    refFunc = composeRefFuncSelect(pref, f);
else
    refFunc = composeRefFuncSelect(pref, f, g);
end
pref.funcheb2.refinementFunction = refFunc;

% Call the funcheb2 constructor:
f = funcheb2(op, vscale, pref);

% Throw a warning:
if ( ~f.epslevel )
    warning('CHEBFUN:FUNCHEB2:compose:convfail',...
        ['Composition with ', func2str(op), ...
        ' failed to converge with ', int2str(length(f)), ' points.']);
end

end

function refFunc = composeRefFuncSelect(pref, f, g)
%COMPOSEREFFUNCSELECT  Choose a refinement function
%   COMPOSEREFFUNCSELECT(PREF, F) chooses a refinement function (either single
%   or double sampling) for constructing OP(F).
%
%   COMPOSEREFFUNCSELECT(PREF, F, G) chooses a refinement function (either
%   single or double sampling) for constructing OP(F, G).

% Double sampling:
if ( strcmp( pref.funcheb2.refinementFunction, 'resampling') )
    if ( nargin == 2 )
        % OP(G1)
        refFunc = @(op, values, pref) composeRefResample1(op, values, pref, f);
    else
        % OP(G1, G2)
        refFunc = @(op, values, pref) composeRefResample2(op, values, pref, f, g);
    end
    
% Single sampling:
else % if ( strcmp( pref.funcheb2.refinementFunction, 'default') )
    if ( nargin == 2 )
        % OP(G1)
        refFunc = @(op, values, pref) composeRefNested1(op, values, pref, f); 
    else
        % OP(G1, G2)
        refFunc = @(op, values, pref) composeRefNested2(op, values, pref, f, g);
    end
end

end

function [values, giveUp] = composeRefResample1(op, values, pref, f)
%COMPOSEREFRESAMPLE1 Refinement function for composing a FUNCHEB2 with resampling.
    
    if ( isempty(values) )
        % Choose initial n based upon minSamples.
        n = 2^ceil(log2(pref.funcheb2.minSamples - 1)) + 1;
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
    if ( n > pref.funcheb2.maxSamples )
        giveUp = true;
        return
    else
        giveUp = false;
    end
    
    % Update f values:
    f = prolong(f, n);
    v1 = f.values;

    % Evaluate the operator
    if ( pref.funcheb2.extrapolate )
        % Avoid evaluating the endpoints:
        valuesTemp = feval(op, v1(2:n-1,:));
        nans = NaN(1, size(valuesTemp, 2));
        values = [nans ; valuesTemp ; nans];
    else
        values = feval(op, v1);
    end
end

function [values, giveUp] = composeRefResample2(op, values, pref, f, g)
%COMPOSEREFRESAMPLE2  Refinement function for composing FUNCHEB2s with resampling.
    
    if ( isempty(values) )
        % Choose initial n based upon minSamples.
        n = 2^ceil(log2(pref.funcheb2.minSamples - 1)) + 1;
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
    if ( n > pref.funcheb2.maxSamples )
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

    if ( pref.funcheb2.extrapolate )
        % Avoid evaluating the endpoints:
        valuesTemp = feval(op, v1(2:n-1,:), v2(2:n-1));
        nans = NaN(1, size(valuesTemp, 2));
        values = [nans ; valuesTemp ; nans];
    else
        values = feval(op, v1, v2);
    end
end

function [values, giveUp] = composeRefNested1(op, values, pref, f)
%COMPOSEREFNESTED1 Refinement function for composing a FUNCHEB2 without resampling.

    if ( isempty(values) )  % We're just starting out:
        [values, giveUp] = composeRefResample1(op, values, pref, f);
        
    else                    % We already have some values
    
        % Compute new n by doubling (we must do this when not resampling).
        n = 2*size(values, 1) - 1;
        
        % n is too large!
        if ( n > pref.funcheb2.maxSamples )
            giveUp = true;
            return
        else
            giveUp = false;
        end
        
        % Update f values:
        f = prolong(f, n);
        v1 = f.values(2:2:end-1,:);

        % Shift the stored values:
        values(1:2:n,:) = values;
        % Compute and insert new ones:
        values(2:2:end-1,:) = feval(op, v1);

    end
end

function [values, giveUp] = composeRefNested2(op, values, pref, f, g)
%COMPOSEREFNESTED2 Refinement function for composing FUNCHEB2s without resampling.
    
    if ( isempty(values) )  % We're just starting out:
        [values, giveUp] = composeRefResample2(op, values, pref, f, g);
        
    else                    % We already have some values
    
        % Compute new n by doubling (we must do this when not resampling).
        n = 2*size(values, 1) - 1;
        
        % n is too large!
        if ( n > pref.funcheb2.maxSamples )
            giveUp = true;
            return
        else
            giveUp = false;
        end
        
        % Update f and g values:
        f = prolong(f, n);
        v1 = f.values(2:2:end-1,:);
        g = prolong(g, n);
        v2 = g.values(2:2:end-1,:);
        
        % Shift the stored values:
        values(1:2:n,:) = values;
        % Compute and insert new ones:
        values(2:2:end-1,:) = feval(op, v1, v2);

    end
end
