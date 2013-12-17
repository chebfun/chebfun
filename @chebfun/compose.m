function f = compose(f, op, g, pref)
%COMPOSE  Composition of CHEBFUN objects.
%   COMPOSE(F, OP) returns a CHEBFUN representing OP(F), where F is also a
%   CHEBFUN object and OP is a function handle.
%
%   COMPOSE(F, OP, G) returns OP(F, G), where F and G are CHEBFUN objects and OP
%   is a function handle. The domains and imensions of F and G should be
%   compatible.
%
%   COMPOSE(F, G) returns a CHEBFUN representing G(F), where both F and G are
%   also CHEBFUN objects. If the range of F is not contained in the domain of G,
%   or if F and G do not have the same dimensions, then an error is thrown.
%
%   COMPOSE(F, OP, PREF), COMPOSE(F, OP, G, PREF), and COMPOSE(F, G, PREF) use
%   the options passed by the CHEBPREF object PREF.
%
%   Note 1: If the locations of required breakpoints in the output are known in
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
% (such as SIN(), for example).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [TODO]: vscale and tolerance?

% Parse inputs:
opIsBinary = false;
if ( (nargin == 4) && ~isempty(g) )           % compose(f, op, g, pref)
    opIsBinary = true;
end
if ( (nargin < 4) || ((nargin == 4) && isempty(pref)) )
    pref = chebpref();
end
if ( nargin == 3 )
    if ( isstruct(g) || isa(g, 'chebpref') )  % compose(f, op, pref)
        pref = chebpref(g);
        g = [];
    else                                      % compose(f, op, g)
        opIsBinary = true;
    end
end
if ( nargin < 3 )                             % compose(f, op) or compose(f, g)
    g = [];
end

%% Special cases:

% There is nothing to do for an empty chebfun!
if ( isempty(f) )
    return
end

if ( isa(op, 'chebfun') )
    % Call the COMPOSETWOCHEBFUNS method if OP is a CHEBFUN object:
    g = op;
    
    if ( numColumns(f) ~= numColumns(g) )
            error('CHEBFUN:compose:dims', 'Matrix dimensions must agree.')
    end
    
    if ( numel(f) == 1 && numel(op) == 1 )
        % Array-valued CHEBFUN case:
        f = composeTwoChebfuns(f, op, pref);
    else
        % QUASIMATRIX case:
        f = cheb2cell(f);
        g = cheb2cell(g);
        for k = numel(f):-1:1
            h(k) = composeTwoChebfuns(f{k}, g{k}, pref);
        end
        f = h;
    end
        
elseif ( opIsBinary )
    % Binary composition:
    
    if ( numColumns(f) ~= numColumns(g) )
            error('CHEBFUN:compose:dims', 'Matrix dimensions must agree.')
    end
    
    if ( numel(f) == 1 && numel(g) == 1 )
        % Array-valued CHEBFUN case:
        f = columnCompose(f, op, g, pref, opIsBinary);
    else
        % QUASIMATRIX case:
        f = cheb2cell(f);
        g = cheb2cell(g);
        for k = numel(f):-1:1
            h(k) = columnCompose(f{k}, op, g{k}, pref, opIsBinary);
        end
        f = h;
    end
    
else
    % Unary composition:
    
    if ( numel(f) == 1 )
        % Array-valued CHEBFUN case:
        f = columnCompose(f, op, g, pref, opIsBinary);
    else
        % QUASIMATRIX case:
        f = cheb2cell(f);
        for k = numel(f):-1:1
            h(k) = columnCompose(f{k}, op, g, pref, opIsBinary);
        end
        f = h;
    end
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = columnCompose(f, op, g, pref, opIsBinary)

%% Initialise:

% Initialise impulses. Note that higher-order impulse data is destroyed by
% compose, so we only require the first impulse 'sheet' (i.e., impulses(:,:,1)):
if ( opIsBinary )
    % Call OVERLAP() if we are composing two CHEBFUN inputs with a binary op:
    [f, g] = overlap(f, g);
    fevalImps = feval(op, f.impulses(:,:,1), g.impulses(:,:,1));
    newImps = fevalImps(1,:);
else
    fevalImps = feval(op, f.impulses(:,:,1));
    newImps = fevalImps(1,:);
end

% Number of piecewise intervals in f:
numInts = numel(f.domain) - 1;

% Initialise storage for the output FUN cell:
newFuns = {};

% Initialise new domain vector:
newDom = f.domain(1);

% Suppress growing vector Mlint warnings (which are inevitable here):
%#ok<*AGROW>

%% Loop through each interval:
for k = 1:numInts

    % Attempt to generate FUN objects using FUN/COMPOSE().
    if ( isempty(g) )
        newFun = compose(f.funs{k}, op, [], pref);
    else
        newFun = compose(f.funs{k}, op, g.funs{k}, pref);
    end
    isHappy = get(newFun, 'ishappy');

    if ( isHappy || ~pref.enableBreakpointDetection )
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

    elseif ( pref.enableBreakpointDetection )

        % If not happy and splitting is on, get a CHEBFUN for that subinterval:
        domk = f.domain(k:k+1);
        if ( opIsBinary )
            newChebfun = chebfun(@(x) feval(op, feval(f, x), feval(g, x)), ...
                domk, pref);
        else
            newChebfun = chebfun(@(x) feval(op, feval(f, x)), domk, pref);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = composeTwoChebfuns(f, g, pref)
%COMPOSETWOCHEBFUNS   Composition of two CHEBFUN objects.
%   COMPOSETWOCHEBFUNS(F, G, PREF) returns the composition of the CHEBFUN
%   objects F and G, G(F) using the CHEBFUN preferences contained in the
%   preference structure PREF.  The range of F must be in the domain of G or
%   else an error is thrown.  An equivalent syntax is G(F).

% Deal with the trivial empty case:
if ( isempty(f) || isempty(g) )
    h = chebfun();
    return
end

%% ERROR CHECKING:

if ( xor(f.isTransposed, g.isTransposed) )
    error('CHEBFUN:composeChebfuns:trans', ...
        'Cannot compose a row CHEBFUN with a column CHEBFUN.');
end

isTransposed = f.isTransposed;
if ( isTransposed )
    % Make everything a column CHEBFUN for now:
    f = transpose(f);
    g = transpose(g);
end

if ( (size(f, 2) > 1) && (size(g, 2) > 1) )
     error('CHEBFUN:composeChebfuns:trans', ...
        'Cannot compose two array-valued CHEBFUN objects.');
end

% f must be a real-valued function:
if ( ~isreal(f) )
    error('CHEBFUN:compose:complex', 'F must be real valued to construct G(F).')
    % warning('CHEBFUN:compose:complex', 'F SHOULD be real valued to construct G(F).');
end

% [TODO]: Requires MINANDMAX().
% % Get epslevels and set a tolerance:
% tol = 10*max(vscale(f).*epslevel(f), vscale(g).*epslevel(g));
% hsf = hscale(f); 
% % Find the range of F:
% mmF = minandmax(f);
% minF = min(mmF(:));
% maxF = max(mmF(:));
% % Range of f must be in the domain of g:
% if ( g.domain(1) > minF + tol*hsf || g.domain(end) < maxF - tol*hsf )
%     error('CHEBFUN:compose:domain', ...
%         'Range of F, [%g, %g], must be in the domain of G, [%g, %g].', ...
%         minF, maxF, g.domain(1), g.domain(end))
% end

% Delta functions:
if ( (size(f.impulses, 3) > 1) || (size(g.impulses, 3) > 1) )
    warning('CHEBFUN:compose:imps',  ...
        'Composition does not handle impulses / delta functions.')
end

%% Locate breakpoints in G:

% If g has breakpoints, find the corresponding x-points in the domain of f:
newDom = f.domain;
if ( numel(g.domain) > 2 )
    gDom = g.domain(2:end-1);
    for k = 1:length(gDom)
        % [TODO]: This requires @CHEBFUN/MINUS.
        % r = roots(f - gDom(k));
        % newDom = [newDom, r(:).']; %#ok<AGROW>
    end
end
newDom = unique(sort(newDom));

% Restrict f to the new domain:
f = restrict(f, newDom);

%% Call COMPOSE():

% Call compose:
h = compose(f, @(f) feval(g, f), pref);

% Fix impulse values:
h.impulses(:,:,1) = feval(g, feval(f, h.domain.'));

if ( isTransposed )
    h = transpose(h);
end

end
