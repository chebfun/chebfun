function f = compose(f, op, g, pref)
%COMPOSE  Compostition of CHEBFUN objects.
%   COMPOSE(F, OP) returns a CHEBFUN representing OP(F) where F is also a
%   CHEBFUN object, and OP is a function handle.
%
%   COMPOSE(F, OP, G) returns OP(F, G) where F and G are CHEBFUN objects, and OP
%   is a function handle.
%
%   COMPOSE(F, G) returns a CHEBFUN representing G(F), where both F and G are
%   also CHEBFUN objects. If the range of F is not in [-1, 1] then an error is
%   thrown.
%
%   COMPOSE(F, OP, PREF), COMPOSE(F, OP, G, PREF), or COMPOSE(F, G, PREF) uses
%   the options passed by the prefences structure PREF.
%
%   Note: If the location of required breakpoints in the output are known in
%   advance, they should be applied to F and/or G using RESTRICT before the call
%   to COMPOSE().
%
% See also RESTRICT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here is a small flowchart of how this process works. (In Chebfun V4 there was
% an additional branch if BLOWUP was ON, but now this is included in
% FUN/COMPOSE.)
%
%                 ----->[compose at FUN level]<-------<-
%                 |     pass/            \fail         |
%                 ^--[next piece]      [splitting?]    |
%                 |                   yes/    \no      |
%                 |   [CHEBFUN constructor]    [fail]  |
%                 |   pass/        \fail         /     |
%                 --[next piece]    \--->[warning]-----^
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [TODO]: vscale?

% Parse inputs:
numChebfuns = 2;
if ( (nargin > 2) && isstruct(g) )
    pref = chebfun.pref(g);
    g = [];
    numChebfuns = 1;
elseif ( (nargin < 4) )
    pref = chebfun.pref();
end
if ( (nargin < 3) || isempty(g) )
    numChebfuns = 1;
    g = [];
end

% There is nothing to do for an empty chebfun!
if ( isempty(f) )
    return
end

% Initialise impulses: (and overlap if there are two chebfuns)
fDom = f.domain(1);
if ( numChebfuns == 1 )
    fImps = feval(op, f.impulses(1,1));
else
    [f, g] = overlap(f, g);
    fImps = feval(op, f.impulses(1,1), g.impulses(1,1));
end

% Number of piecewise smooth components:
numFuns = numel(f.domain) - 1;

% Initialise storage for output funs:
fFuns = {};

% Merge preferences with fun.prefs:
pref = fun.pref(pref, pref.chebfun);

% Loop through the funs:
for k = 1:numFuns

    % Attempt to generate FUN objects using the FUN/COMPOSE().
    if ( isempty(g) )
        newFun = compose(f.funs{k}, op, [], pref);
    else
        newFun = compose(f.funs{k}, op, g.funs{k}, pref);
    end
    isHappy = get(newFun, 'ishappy');

    if ( isHappy || ~pref.chebfun.splitting )
        % If we're happy or not allowed to split, this will do.

        if ( ~isHappy )
            % Throw a warning if we're not happy:
            try
                str = ['with function ' func2str(op)];
            catch ME %#ok<NASGU>
                str = '';
            end
            warning('CHEBFUN:compose:resolve', ['Composition ', str, ...
                ' not resolved using ' int2str(length(newFun)), ...
                ' points. Have you tried ''splitting on''?']);
        end

        % Store new fun:
        fFuns = [fFuns, {newFun}];                    %#ok<AGROW>
        % Store new ends:
        fDom = [fDom, f.domain(k+1)];                 %#ok<AGROW>
        % Store new imps:
        if ( isempty(g) )
            fImps = [fImps, feval(op, f.impulses(1,k+1))]; %#ok<AGROW>
        else
            fImps = [fImps, feval(op, f.impulses(1,k+1), g.impulses(1,k+1))]; %#ok<AGROW>
        end

    elseif ( pref.chebfun.splitting )

        % If not happy and splitting is on, get a chebfun for that subinterval.
        domk = f.domain(k:k+1);
        if ( numChebfuns == 1 )
            newCfun = chebfun(@(x) feval(op, feval(f, x)), domk, pref);
        else
            newCfun = chebfun(@(x) feval(op, feval(f, x), feval(g, x)), ...
                domk, pref);
        end

        isHappy = get(newCfun, 'ishappy');
        if ( ~isHappy )
            % Throw a warning if we're not happy:
            try
                str = ['with function ' func2str(op)];
            catch ME %#ok<NASGU>
                str = '';
            end
            warning('CHEBFUN:compose:resolve', ['Composition ', str, ...
                ' not resolved using ' int2str(length(newFun)) ' points.']);
        end

        % Store new funs:
        fFuns = [fFuns, newCfun.funs];                %#ok<AGROW>
        % Store new ends:
        fDom = [fDom, newCfun.domain(2:end)];     %#ok<AGROW>
        % Store new imps:
        fImps = [fImps, newCfun.impulses(1,2:end)];     %#ok<AGROW>

    end

end

% Put the funs back into a chebfun:
f.funs = fFuns;
f.domain = fDom;
f.impulses = fImps;

end
