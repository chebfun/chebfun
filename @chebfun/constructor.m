function [funs, ends] = constructor(op, domain, pref)
%CONSTRUCTOR   CHEBFUN constructor.
%   FUNS = CONSTRUCTOR(OP, DOMAIN) constructs the piecewise components (known as
%   "FUNS") used by a CHEBFUN object to represent the function OP on the
%   interval DOMAIN. OP must be a function_handle, string, numerical vector, or
%   a cell array containing a combination of these first three data types. In
%   the later case, the number of elements in the array must be one less than
%   the length of the DOMAIN vector.
%
%   If OP is a function_handle or a string, it should be vectorised in that it
%   accepts a column vector of length N and return a matrix of size N x M. If M
%   ~= 1, we say the resulting CHEBFUN is "array-valued".
%
%   CONSTRUCTOR(OP, DOMAIN, PREF), where PREF is a structure returned by 
%   CHEBFUN.PREF(), allows alternative construction preferences to be passed to
%   the constructor. See >> help chebfun/pref for more details on preferences.
%
%   In particular, if PREF.CHEBFUN.SPLITTING = TRUE and OP is a function_handle
%   or a string, then the constructor adaptively introduces additional
%   breakpoints into the domain so as to better represent the function. These
%   are returned as the second output argument in [FUNS, NEWDOMAIN] =
%   CONSTRUCTOR(OP, DOMAIN).
%
% See also CHEBFUN, CHEBFUN/PREF.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Initial setup:
ends = domain;
numIntervals = numel(ends) - 1;

% Initialise hscale and vscale:
hscale = norm(domain, inf);
if ( isinf(hscale) )
    hscale = 1;
end
vscale = 0;

% Sanity check:
if ( iscell(op) && (numel(op) ~= numIntervals) )
    error('CHEBFUN:constructor:cellInput', ...
        ['Number of cell elements in OP must match the number of ', ...
         'intervals in DOMAIN.'])
end    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%  SPLITTING OFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In 'OFF' mode, seek only one piece with length < maxdegree.
if ( ~pref.chebfun.splitting )
    % Set maximum number of sample points:
    maxn = pref.chebfun.maxdegree + 1;
    % Initialise the FUN array:
    funs{numIntervals} = fun.constructor();
    % Merge with FUN preferences:
    pref = fun.pref(pref, pref.chebfun);
    % We only want to throw this warning once:
    warningThrown = false;
    % Loop over the intervals:
    for k = 1:numIntervals
        endsk = ends(k:k+1);
        % Unwrap if OP is a cell:
        if ( iscell(op) )
            opk = op{k};
        else
            opk = op;
        end
        % Call GETFUN() (which calls FUN.CONSTRUCTOR()):
        [funs{k}, ishappy, vscale] = getFun(opk, endsk, vscale, hscale, pref);
        % Warn if unhappy (as we're unable to split the domain to improve):
        if ( ~ishappy && ~warningThrown)
            warning('CHEBFUN:constructor', ...
                ['Function not resolved using %d pts.', ...
                ' Have you tried ''splitting on''?'], maxn);
            warningThrown = true;
        end
    end
    return
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%  SPLITTING ON %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In 'ON' mode, seek only many pieces with total length < maxlength.

% Set the maximum degree:
pref.chebfun.maxSamples = pref.chebfun.splitdegree;
% We extrapolate when splitting so that we can construct functions like
% chebfun(@sign,[-1 1]), which otherwise would not be happy at x = 0.
pref.chebfun.extrapolate = true;
% Merge with FUN preferences:
pref = fun.pref(pref, pref.chebfun);

% Initialise the FUN array:
funs{numIntervals} = fun.constructor();
% Initialise happiness:
ishappy = ones(1, numel(ends) - 1);

% Try to get one smooth piece for the entire domain before splitting:
for k = 1:numIntervals
    % Unwrap if OP is a cell:
    if ( iscell(op) )
        opk = op{k};
    else
        opk = op;
    end
    [funs{k}, ishappy(k), vscale] = ...
        getFun(opk, ends(k:k+1), vscale, hscale, pref);
end
sad = ~ishappy;

% MAIN LOOP. If the above didn't work, enter main loop and start splitting.
% (Note that at least one new breakpoint will be introduced).
while ( any(sad) )
    % If a FUN is sad in a subinterval, split this subinterval.

    % Choose a subinterval to improve:
%     % Old choice = the first sad interval:
%     k = find(sad, 1, 'first');
    % New choice = the largest sad interval:
    diffEnds = diff(ends);
    diffEnds(~sad) = 0;
    [ignored, k] = max(diffEnds);
    
    % Ends of this subinterval:
    a = ends(k);
    b = ends(k+1);
    
    % Unwrap if OP is a cell:
    if ( iscell(op) )
        opk = op{k};
    else
        opk = op;
    end

    % Locate an edge/split location:
    edge = chebfun.detectEdge(opk, [a, b], vscale, hscale);

    % Try to obtain happy child FUN objects on each new subinterval:
    [childLeft, happyLeft, vscale] = ...
        getFun(opk, [a, edge], vscale, hscale, pref);
    [childRight, happyRight, vscale] = ...
        getFun(opk, [edge, b], vscale, hscale, pref);

    % Check for happiness/sadness:
    sad = [sad(1:k-1), ~happyLeft, ~happyRight, sad(k+1:end)];

    % Insert new pieces in to existing funs:
    funs = [funs(1:k-1), {childLeft, childRight}, funs(k+1:end)];
    ends = [ends(1:k), edge, ends(k+1:end)];

    % If a cell was given, we must store pieces on new intervals:
    if ( iscell(op) )
        op = [op(1:k), {opk}, op(k+1:end)];
    end

    % Fail if too many points are required:
    if ( sum(cellfun(@length, funs)) > pref.chebfun.maxlength )
        warning('Function not resolved using %d pts.', ...
            sum(cellfun(@length, funs)));
        return
    end

end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GETFUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [g, ishappy, vscale] = getFun(op, domain, vscale, hscale, pref)
%GETFUN controls the construction of funs

% If the interval is very small then skip adaptation and treat OP as a constant:
if ( diff(domain) < 2e-14*hscale )
    op = op(mean(domain));
end

% Call the fun constructor:
g = fun.constructor(op, domain, vscale, hscale, pref);
% See if the construction was happy:
ishappy = get(g, 'ishappy');
% Update the vertical scale:
vscale = max(vscale, get(g, 'vscale'));

end
