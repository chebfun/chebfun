function f = compose(f, op, g, pref)
%COMPOSE  Compostition of CHEBTECH2 objects.
%   COMPOSE(F, OP) returns a CHEBTECH2 representing OP(F) where F is also a
%   CHEBTECH2 object, and OP is a function handle.
%
%   COMPOSE(F, OP, G) returns OP(F, G) where F and G are CHEBTECH2 objects,
%   and OP is a function handle.
%
%   COMPOSE(F, G) returns a CHEBTECH2 representing G(F), where both F and G are
%   also CHEBTECH objects. If the range of F is not in [-1, 1] then an error is
%   thrown.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This file defines the refinement functions which are called by the constructor
% at the @chebtech level. There are refinement functions for both nested and
% resampled grids, both for compositions of the form op(f) and op(f, g).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nfuns = 2;
% Parse inputs:
if ( nargin > 2 && isstruct(g) )
    pref = g;
    g = [];
    nfuns = 1;
elseif ( nargin < 4 )
    pref = chebtech.pref();
end
if ( nargin < 3 || isempty(g) )
    nfuns = 1;
    g = [];
end

% Choose a sampling strategy:
if ( ~ischar(pref.chebtech.refinementFunction) )
    % A user-defined refinement has been passed.
    refFunc = pref.chebtech.refinementFunction;
elseif ( strcmp( pref.chebtech.refinementFunction, 'resampling') )
    if ( nfuns == 1 )   % OP(G1) resampling.
        refFunc = @(op, values, pref) composeResample1(op, values, pref, f);
    else                % OP(G1, G2) resampling.
        refFunc = @(op, values, pref) composeResample2(op, values, pref, f, g);
    end    
else
    if ( nfuns == 1 )   % OP(G1) nested.
        refFunc = @(op, values, pref) composeNested1(op, values, pref, f); 
    else                % OP(G1, G2) nested.
        refFunc = @(op, values, pref) composeNested2(op, values, pref, f, g);
    end
end
% Assign to preference structure:
pref.chebtech.refinementFunction = refFunc;

% Call parent compose:
f = compose@chebtech(f, op, g, pref);

end

function [values, giveUp] = composeResample1(op, values, pref, f)
%COMPOSERESAMPLE1 Refinement function for composing a CHEBTECH2 with resampling.
    
    if ( isempty(values) )
        % Choose initial n based upon minSamples.
        n = 2^ceil(log2(pref.chebtech.minSamples - 1)) + 1;
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
    if ( n > pref.chebtech.maxSamples )
        giveUp = true;
        return
    else
        giveUp = false;
    end
    
    % Update f values:
    f = prolong(f, n);
    v1 = f.values;

    % Evaluate the operator
    if ( pref.chebtech.extrapolate )
        % Avoid evaluating the endpoints:
        valuesTemp = feval(op, v1(2:n-1,:));
        nans = NaN(1, size(valuesTemp, 2));
        values = [nans ; valuesTemp ; nans];
    else
        values = feval(op, v1);
    end
end

function [values, giveUp] = composeResample2(op, values, pref, f, g)
%COMPOSERESAMPLE2  Refinement function for composing CHEBTECH2s with resampling.
    
    if ( isempty(values) )
        % Choose initial n based upon minSamples.
        n = 2^ceil(log2(pref.chebtech.minSamples - 1)) + 1;
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
    if ( n > pref.chebtech.maxSamples )
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

    if ( pref.chebtech.extrapolate )
        % Avoid evaluating the endpoints:
        valuesTemp = feval(op, v1(2:n-1,:), v2(2:n-1));
        nans = NaN(1, size(valuesTemp, 2));
        values = [nans ; valuesTemp ; nans];
    else
        values = feval(op, v1, v2);
    end
end

function [values, giveUp] = composeNested1(op, values, pref, f)
%COMPOSENESTED1 Refinement function for composing a CHEBTECH2 without resampling.

    if ( isempty(values) )  % We're just starting out:
        [values, giveUp] = composeResample1(op, values, pref, f);
        
    else                    % We already have some values
    
        % Compute new n by doubling (we must do this when not resampling).
        n = 2*size(values, 1) - 1;
        
        % n is too large!
        if ( n > pref.chebtech.maxSamples )
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

function [values, giveUp] = composeNested2(op, values, pref, f, g)
%COMPOSENESTED2 Refinement function for composing CHEBTECH2s without resampling.
    
    if ( isempty(values) )  % We're just starting out:
        [values, giveUp] = composeResample2(op, values, pref, f, g);
        
    else                    % We already have some values
    
        % Compute new n by doubling (we must do this when not resampling).
        n = 2*size(values, 1) - 1;
        
        % n is too large!
        if ( n > pref.chebtech.maxSamples )
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
