function f = compose(f, op, g, pref)
%COMPOSE   Composition of CHEBTECH1 objects.
%   COMPOSE(F, OP) returns a CHEBTECH1 representing OP(F), where F is also a
%   CHEBTECH1 object, and OP is a function handle.
%
%   COMPOSE(F, OP, G) returns a CHEBTECH1 representing OP(F, G), where F and G
%   are CHEBTECH objects, and OP is a function handle.
%
%   COMPOSE(F, G) returns a CHEBTECH1 representing G(F), where both F and G are
%   also CHEBTECH objects. If the range of F is not in [-1, 1] then an error is
%   thrown.
%
%   COMPOSE(F, OP, G, PREF) or COMPOSE(F, OP, [], PREF) uses the options passed
%   by the preferences structure PREF to build the returned CHEBTECH1.  In
%   particular, one can set PREF.REFINEMENTFUNCTION to be a function which takes
%   advantage of F and possibly OP or G being CHEBTECH objects.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This file defines the refinement functions which are called by the
% constructor at the @chebtech level. There are refinement functions for
% resampled grids only, both for compositions of the form op(f) and op(f, g).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse inputs:
if ( nargin < 4 )
    pref = f.techPref();
else
    pref = f.techPref(pref);
end

if ( (nargin < 3) || isempty(g) )
    nfuns = 1;
    g = [];
else
    nfuns = 2;
end

% Choose a sampling strategy:
if ( ~ischar(pref.refinementFunction) )
    % A user-defined refinement has been passed.
    refFunc = pref.refinementFunction;
else
    if ( nfuns == 1 )   % OP(G1) resampling.
        refFunc = @(op, values, pref) composeResample1(op, values, pref, f);
    else                % OP(G1, G2) resampling.
        refFunc = @(op, values, pref) composeResample2(op, values, pref, f, g);
    end    
end

% Assign to preference structure:
pref.refinementFunction = refFunc;

% Call superclass COMPOSE:
f = compose@chebtech(f, op, g, pref);

end

function [values, giveUp] = composeResample1(op, values, pref, f)
%COMPOSERESAMPLE1   Refinement function for composing a CHEBTECH1 with
%resampling.
    
    if ( isempty(values) )
        % Choose initial n based upon minSamples.
        n = 2^ceil(log2(pref.minSamples - 1)) + 1;
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
    if ( n > pref.maxSamples )
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
%COMPOSERESAMPLE2   Refinement function for composing CHEBTECH1 objects with
%resampling.
    
    if ( isempty(values) )
        % Choose initial n based upon minSamples.
        n = 2^ceil(log2(pref.minSamples - 1)) + 1;
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
    if ( n > pref.maxSamples )
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
