function f = compose(f, op, g, pref)
%COMPOSE  Compostition of CHEBFUN objects.
%   COMPOSE(F, OP) returns a CHEBFUN representing OP(F), where F is also a
%   CHEBFUN object and OP is a function handle.
%
%   COMPOSE(F, OP, G) returns OP(F, G), where F and G are CHEBFUN objects and OP
%   is a function handle. The domains of F and G should be compatible.
%
%   COMPOSE(F, G) returns a CHEBFUN representing G(F), where both F and G are
%   also CHEBFUN objects. If the range of F is not in [-1, 1] then an error is
%   thrown.
%
%   COMPOSE(F, OP, PREF), COMPOSE(F, OP, G, PREF), and COMPOSE(F, G, PREF) use
%   the options passed by the prefences structure PREF.
%
%   Note 1: If the location of required breakpoints in the output are known in
%   advance, they should be applied to F and/or G using RESTRICT() before the
%   call to COMPOSE().
%
%   Note 2: Any higher-order impulse/delta function data in F and/or G is
%   ignored by COMPOSE().

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
% It is not expected that the user calls COMPOSE() directly. It is usually
% accessed via other @CHEBFUN methods which use it for their implementation
% (such a SIN()).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [TODO]: vscale and tolerance?

%% Parse inputs:
numChebfuns = 2;
if ( (nargin == 3) && isstruct(g) ) % compose(f, op, pref)
    pref = chebfun.pref(g);
    g = [];
    numChebfuns = 1;
elseif ( (nargin < 4) )             % compose(f, op, g)
    pref = chebfun.pref();
end
if ( (nargin < 3) || isempty(g) )   % compose(f, op) or compose(f, op, [], pref)
    numChebfuns = 1;
    g = [];
end

%% Special cases:

% There is nothing to do for an empty chebfun!
if ( isempty(f) )
    return
end

% Call the COMPOSECHEBFUNS method if OP is a CHEBFUN object:
if ( isa(op, 'chebfun') )
    f = composeChebfuns(f, g, pref);
    return
end

%% Initialise:

% Initialise impulses. Note that higher-order impulse data is destroyed by
% compose, so we only require the first impulse 'sheet' (i.e., impulses(:,:,1)):
if ( numChebfuns == 1 )
    fevalImps = feval(op, f.impulses(:,:,1));
    newImps = fevalImps(1,:);
else
    % Call OVERLAP() if there are two CHEBFUN inputs:
    [f, g] = overlap(f, g);
    fevalImps = feval(op, f.impulses(:,:,1), g.impulses(:,:,1));
    newImps = fevalImps(1,:);
end

% Number of piecewise intervals in f:
numInts = numel(f.domain) - 1;

% Initialise storage for the output FUN cell:
newFuns = {};

% Initialise new domain vector:
newDom = f.domain(1);

% Merge preferences with FUN.PREF();
pref = fun.pref(pref, pref.chebfun);

% Suppress growing vector Mlint warnings (which are inevitable here):
%#ok<*AGROW>

%% Loop through the each interval:
for k = 1:numInts

    % Attempt to generate FUN objects using FUN/COMPOSE().
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
                str = ['with function ', func2str(op)];
            catch ME %#ok<NASGU>
                str = '';
            end
            warning('CHEBFUN:compose:resolve', ['Composition ', str, ...
                ' not resolved using ', int2str(length(newFun)), ...
                ' points. Have you tried ''splitting on''?']);
        end

        % Store new FUN in cell array:
        newFuns = [newFuns, {newFun}];                   
        % Store new ends:
        newDom = [newDom, f.domain(k+1)];
        % Store new impulses: (Note, will only be a matrix - not a tensor)
        if ( isempty(g) )
            newImps = [newImps ; fevalImps(k+1,:)];
        else
            newImps = [newImps ; fevalImps(k+1,:)];
        end

    elseif ( pref.chebfun.splitting )

        % If not happy and splitting is on, get a CHEBFUN for that subinterval:
        domk = f.domain(k:k+1);
        if ( numChebfuns == 1 )
            newChebfun = chebfun(@(x) feval(op, feval(f, x)), domk, pref);
        else
            newChebfun = chebfun(@(x) feval(op, feval(f, x), feval(g, x)), ...
                domk, pref);
        end

        if ( ~get(newChebfun, 'ishappy') ) % Throw a warning if we're not happy:
            try
                str = ['with function ', func2str(op)];
            catch ME %#ok<NASGU>
                str = '';
            end
            warning('CHEBFUN:compose:resolve', ['Composition ', str, ...
                ' not resolved using ', int2str(length(newChebfun)), ...
                ', points.']);
        end

        % Store new FUN objects:
        newFuns = [newFuns, newChebfun.funs];
        % Store new ends:
        newDom = [newDom, newChebfun.domain(2:end)];
        % Store new impulses; (Note, will only be a matrix - not a tensor)
        newImps = [newImps ; newChebfun.impulses(2:end-1,:) ; ...
            fevalImps(k+1,:) ];

    end

end

%% Prepare output:

% Put the FUN cell, domain, and impulses back into a CHEBFUN:
f.funs = newFuns;
f.domain = newDom;
f.impulses = newImps;

end
