function f = compose(f, op, g, pref)
%COMPOSE   Composition of CHEBTECH1 objects.
%   COMPOSE(F, OP) returns a CHEBTECH1 representing OP(F) where F is also a
%   CHEBTECH1 object, and OP is a function handle.
%
%   COMPOSE(F, OP, G) returns OP(F, G) where F and G are CHEBTECH1 objects,
%   and OP is a function handle.
%
%   COMPOSE(F, G) returns a CHEBTECH1 representing G(F), where both F and G are
%   also CHEBTECH objects. If the range of F is not in [-1, 1] then an error is
%   thrown.
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
    pref = chebtech.pref();
end
if ( (nargin < 3) || isempty(g) )
    nfuns = 1;
    g = [];
end

% Choose a sampling strategy:
if ( ~ischar(pref.chebtech.refinementFunction) )
    % A user-defined refinement has been passed.
    refFunc = pref.chebtech.refinementFunction;
else
    if ( nfuns == 1 )   % OP(G1) resampling.
        refFunc = @(op, values, pref) composeResample1(op, values, pref, f);
    else                % OP(G1, G2) resampling.
        refFunc = @(op, values, pref) composeResample2(op, values, pref, f, g);
    end    
end

% Assign to preference structure:
pref.chebtech.refinementFunction = refFunc;

% Call superclass COMPOSE:
f = compose@chebtech(f, op, g, pref);

end

function [values, giveUp] = composeResample1(op, values, pref, f)
%COMPOSERESAMPLE1   Refinement function for composing a CHEBTECH1 with resampling.
    
    if ( isempty(values) )
        % Choose initial n based upon minSamples.
        n = 2^ceil(log2(pref.chebtech.minSamples - 1)) + 1;
    else
        % (Approximately) powers of sqrt(2):
        pow = log2(size(values, 1) - 1);
        if ( (pow == floor(pow)) && (pow > 5) )
            n = round(2^(floor(pow) + .5)) + 1;
            n = n - mod(n, 2) + 1;
        else
            n = 2^(floor(pow) + 1) + 1;
        end
    end
    
    % n is too large.
    if ( n > pref.chebtech.maxSamples )
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

function [values, giveUp] = composeResample2(op, values, pref, f, g)
%COMPOSERESAMPLE2   Refinement function for composing CHEBTECH1s with resampling.
    
    if ( isempty(values) )
        % Choose initial n based upon minSamples.
        n = 2^ceil(log2(pref.chebtech.minSamples - 1)) + 1;
    else
        % (Approximately) powers of sqrt(2):
        pow = log2(size(values, 1) - 1);
        if ( (pow == floor(pow)) && (pow > 5) )
            n = round(2^(floor(pow) + .5)) + 1;
            n = n - mod(n, 2) + 1;
        else
            n = 2^(floor(pow) + 1) + 1;
        end
    end
    
    % n is too large:
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
    values = feval(op, v1, v2);

end
