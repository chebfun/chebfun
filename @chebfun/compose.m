function f = compose(f, op, g, pref)
%COMPOSE  Composition of CHEBFUN objects.
%   COMPOSE(F, OP) returns a CHEBFUN representing OP(F), where F is also a
%   CHEBFUN object and OP is a function handle.
%
%   COMPOSE(F, OP, G) returns OP(F, G), where F and G are CHEBFUN objects 
%   and OP is a function handle. The domains and dimensions of F and G 
%   should be compatible.
%
%   COMPOSE(F, G) returns a CHEBFUN representing G(F), where both F and G 
%   are also CHEBFUN objects. If the range of F is not contained in the 
%   domain of G, or if F and G do not have the same dimensions, then an 
%   error is thrown.
%
%   COMPOSE(F, G) where G is a CHEBFUN3 or CHEBFUN3V and F has three real
%   columns returns a CHEBFUN representing G(F).  Similarly if G is a 
%   CHEBFUN2 or CHEBFUN2V and F has two real columns.  If F has one column,
%   then G(F) is iterpreted as G(real(F), imag(F)), regardless of whether F
%   is real or complex.
%
%   COMPOSE(F, OP, PREF), COMPOSE(F, OP, G, PREF), and COMPOSE(F, G, PREF) 
%   use the options passed by the CHEBFUNPREF object PREF.
%
%   Note: If the locations of required breakpoints in the output are known 
%   in advance, they should be applied to F and/or G using RESTRICT() 
%   before the call to COMPOSE().

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

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

% [TODO]: vscale?

% [TODO]: 
%  The lower level constructors and compose() methods all accept a 'data'
%  structure input. In the long term it might be useful to allow the chebfun
%  compose() methods also to accept this (note that chebfun.constructor())
%  already does). However, this is a nontrivial amount of work, and currently
%  there's nothing that would make use of this feature 9for example, a
%  composition whereby you know the form of the exponents in the solution). If
%  that changes, we should revisit this issue.

% Parse inputs:
opIsBinary = false;
if ( (nargin == 4) && ~isempty(g) )           % compose(f, op, g, pref)
    opIsBinary = true;
end
if ( (nargin < 4) || ((nargin == 4) && isempty(pref)) )
    pref = chebfunpref();
end
if ( nargin == 3 )
    if ( isstruct(g) || isa(g, 'chebfunpref') )  % compose(f, op, pref)
        pref = chebfunpref(g);
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
    
    if ( numColumns(f) > 1 && numColumns(g) > 1 )
        error('CHEBFUN:CHEBFUN:compose:trans', ...
            'Cannot compose two array-valued CHEBFUN objects.');
    end
    
    if ( numel(f) == 1 && numel(op) == 1 )
        % Array-valued CHEBFUN case:
        f = composeTwoChebfuns(f, op, pref);
    else
        % QUASIMATRIX case:
        f = cheb2cell(f);
        g = cheb2cell(g);
        if ( numel(f) > 1 )
            for k = numel(f):-1:1
                h(k) = composeTwoChebfuns(f{k}, g{1}, pref);
            end
        else
            for k = numel(g):-1:1
                h(k) = composeTwoChebfuns(f{1}, g{k}, pref);
            end
        end
        f = h;
    end
    
elseif ( isa(op, 'chebfun2') )
    if ( size(f, 2) == 1 )
        % Interpret as OP(real(f), imag(f)), regardless whether f is real- or
        % complex-valued.
        f = compose([ real(f), imag(f) ], op);
    elseif ( size(f, 2) == 2 )
        % Extract the two columns of f:
        x = extractColumns(f, 1);
        y = extractColumns(f, 2);
        
        % Make sure each column is real-valued:
        if ( ~isreal(x) || ~isreal(y) )
            error('CHEBFUN:CHEBFUN:compose:complex2', ...
                'CHEBFUN2(f) is defined for a chebfun f with one complex or two real columns.')
        end
        
        % Compose:
        if ( isPeriodicTech(x) && isPeriodicTech(y) )
            % OP(f) should be periodic if f is:
            f = chebfun(@(t) op(feval(x, t), feval(y, t)), x.domain, 'trig');
        else
            f = chebfun(@(t) op(feval(x, t), feval(y, t)), x.domain);
        end
    else
        error('CHEBFUN:CHEBFUN:compose:Cheb2ofCheb', ...
            'CHEBFUN2(f) is defined for a chebfun f with one complex or two real columns.')
    end
    
elseif ( isa(op, 'chebfun2v') )
    if ( size(f, 2) == 1 )
        % Interpret as OP(real(f), imag(f)), regardless whether f is real- or
        % complex-valued.
        f = compose([ real(f), imag(f) ], op);
    elseif ( size(f, 2) == 2 )
        % Call compose for each component of OP.
        F = compose(f, op(1));
        for jj = 2:op.nComponents
            F = [ F, compose(f, op(jj)) ];
        end
        f = F;
    else
        error('CHEBFUN:CHEBFUN:compose:Cheb2VofCheb', ...
            'CHEBFUN2V(f) is defined for a chebfun f with one complex or two real columns.')
    end
    
elseif ( isa(op, 'chebfun3') )
    if ( size(f, 2) == 3 )
        % Extract the columns of f:
        x = extractColumns(f, 1);
        y = extractColumns(f, 2);
        z = extractColumns(f, 3);
        
        % Make sure all components are real-valued.
        if ( ~isreal(x) || ~isreal(y) || ~isreal(z) )
            error('CHEBFUN:CHEBFUN:compose:complex3', ...
                'CHEBFUN3(F) is defined for a chebfun F with three real columns.')
        end
        
        % Compose:
        if ( isPeriodicTech(x) && isPeriodicTech(y) && isPeriodicTech(z) )
            % OP(f) should be periodic if f is:
            f = chebfun(@(t) op(feval(x, t), feval(y, t), feval(z, t)), ...
                x.domain, 'trig');
        else
            f = chebfun(@(t) op(feval(x, t), feval(y, t), feval(z, t)), ...
                x.domain);
        end
    else
        error('CHEBFUN:CHEBFUN:compose:Cheb3ofCheb', ...
            'CHEBFUN3(F) is defined for a chebfun F with three columns.')
    end
    
elseif ( isa(op, 'chebfun3v') )
    if ( size(f, 2) == 3 )
        % Call compose for each component of OP.
        F = compose(f, op(1));
        for jj = 2:op.nComponents
            F = [ F, compose(f, op(jj)) ];
        end
        f = F;
    else
       error('CHEBFUN:CHEBFUN:compose:Cheb3VofCheb', ...
            'CHEBFUN3V(F) is defined for a chebfun F with three columns.')
    end 
        
elseif ( opIsBinary )
    % Binary composition:
    
    dimCheck(f, g);
    
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

% Initialise pointValues.

if ( opIsBinary )
    % Call OVERLAP() if we are composing two CHEBFUN inputs with a binary op:
    [f, g] = overlap(f, g);
    pointValues = feval(op, f.pointValues, g.pointValues);
    newPointValues = pointValues(1,:);
else
    pointValues = feval(op, f.pointValues);
    newPointValues = pointValues(1,:);
end

% Number of piecewise intervals in f:
numInts = numel(f.domain) - 1;

% Initialise storage for the output FUN cell:
newFuns = {};

% Initialise new domain vector:
newDom = f.domain(1);

if ( pref.splitting) % Is splitting on?
    % Set the maximum length (i.e., number of sample points for CHEBTECH):
    pref.techPrefs.maxLength = pref.splitPrefs.splitLength;
end

if ( numInts > 1 )
    % Force extrapolation if there are multiple pieces:
    pref.extrapolate = 1;
end

% Make sure that the tech preference is set correctly.
%
% TODO:  What is the right thing to do here when composing two chebfuns f(g)?
% This line assumes that the compose method for f is the one that gets called
% ultimately instead of the one for g (which is true currently but could
% change).
pref.tech = get(f.funs{1}, 'tech');

% Suppress growing vector Mlint warnings (which are inevitable here):
%#ok<*AGROW>

%% Loop through each interval:
for k = 1:numInts
    
    % Attempt to generate FUN objects using FUN/COMPOSE().
    if ( isempty(g) )
        newFun = compose(f.funs{k}, op, [], [], pref);
    else
        newFun = compose(f.funs{k}, op, g.funs{k}, [], pref);
    end
    isHappy = get(newFun, 'ishappy');

    if ( isHappy || ~pref.splitting )
        % If we're happy or not allowed to split, this will do.

        if ( ~isHappy )
            % Throw a warning if we're not happy:
            try
                str = ['with function ', func2str(op)];
            catch ME %#ok<NASGU>
                str = '';
            end
            warning('CHEBFUN:CHEBFUN:compose:resolve', ['Composition ', str, ...
                ' not resolved using ', int2str(length(newFun)), ...
                ' points. Have you tried ''splitting on''?']);
        end

        % Store new FUN in cell array:
        newFuns = [newFuns, {newFun}];                   
        % Store new ends:
        newDom = [newDom, f.domain(k+1)];
        % Store new pointValues: (Note, will only be a matrix - not a tensor)
        if ( isempty(g) )
            newPointValues = [newPointValues ; pointValues(k+1,:)];
        else
            newPointValues = [newPointValues ; pointValues(k+1,:)];
        end

    elseif ( pref.splitting )

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
            warning('CHEBFUN:CHEBFUN:compose:resolve', ['Composition ', str, ...
                ' not resolved using ', int2str(length(newChebfun)), ...
                ', points.']);
        end

        % Store new FUN objects:
        newFuns = [newFuns, newChebfun.funs];
        % Store new ends:
        newDom = [newDom, newChebfun.domain(2:end)];
        % Store new pointValues; (Note, will only be a matrix - not a tensor)
        newPointValues = [newPointValues ; newChebfun.pointValues(2:end-1,:) ; ...
            pointValues(k+1,:) ];

    end

end

%% Prepare output:

% Put the FUN cell, domain, and pointValues back into a CHEBFUN:
f.funs = newFuns;
f.domain = newDom;
f.pointValues = newPointValues;

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
    error('CHEBFUN:CHEBFUN:compose:composeTwoChebfuns:trans', ...
        'Cannot compose a row CHEBFUN with a column CHEBFUN.');
end

isTransposed = f.isTransposed;
if ( isTransposed )
    % Make everything a column CHEBFUN for now:
    f = transpose(f);
    g = transpose(g);
end

% f must be a real-valued function:
if ( ~isreal(f) )
%     error('CHEBFUN:CHEBFUN:compose:complex', 'F must be real valued to construct G(F).')
    warning('CHEBFUN:CHEBFUN:compose:composeTwoChebfuns:complex', ...
        ['F should be real valued to construct G(F).\n', ...
         'Results may be inaccurate if G is not a polynomial.']);
     warning('off', 'CHEBFUN:CHEBFUN:compose:composeTwoChebfuns:complex');
else

    % Set a tolerance:
    tol = 1000*eps*max(vscale(f), vscale(g));
    hsf = hscale(f); 
    % Find the range of F:
    mmF = minandmax(f);
    minF = min(mmF(:));
    maxF = max(mmF(:));
    % Range of f must be in the domain of g:
    if ( g.domain(1) > minF + tol*hsf || g.domain(end) < maxF - tol*hsf )
        error('CHEBFUN:CHEBFUN:compose:domain', ...
            'Range of F, [%g, %g], must be in the domain of G, [%g, %g].', ...
            minF, maxF, g.domain(1), g.domain(end))
    end
end

if ( isdelta(f) || isdelta(g) )
    warning('CHEBFUN:CHEBFUN:compose:composeTwoChebfuns:deltas', ...
        'Composition ignores delta functions. Results may not make any sense.');
end

%% Locate breakpoints in G:

% If g has breakpoints, find the corresponding x-points in the domain of f:
if ( numel(g.domain) > 2 )
    newDom = f.domain;
    gDom = g.domain(2:end-1);
    for k = 1:length(gDom)
        r = roots(f - gDom(k));
        newDom = [newDom, r(:).'];
    end
    newDom = chebfun.tolUnique(sort(newDom));
    
    % Restrict f to the new domain:
    f = restrict(f, newDom);
end

%% Call COMPOSE():

% Call compose:
h = compose(f, @(f) feval(g, f), pref);

% Fix impulse values:
h.pointValues = feval(g, feval(f, h.domain.'));

if ( isTransposed )
    h = transpose(h);
end

end
