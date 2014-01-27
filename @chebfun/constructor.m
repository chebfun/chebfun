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
%   CONSTRUCTOR(OP, DOMAIN, PREF), where PREF is a CHEBPREF object, allows
%   alternative construction preferences to be passed to the constructor. See
%   >> help chebpref for more details on preferences.
%
%   In particular, if PREF.ENABLEBREAKPOINTDETECTION = TRUE and OP is a
%   function_handle or a string, then the constructor adaptively introduces
%   additional breakpoints into the domain so as to better represent the
%   function. These are returned as the second output argument in [FUNS, END] =
%   CONSTRUCTOR(OP, DOMAIN).
%
% See also CHEBFUN, CHEBPREF.
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
vscale = pref.scale;

% Sanity check:
if ( iscell(op) && (numel(op) ~= numIntervals) )
    error('CHEBFUN:constructor:cellInput', ...
        ['Number of cell elements in OP must match the number of ', ...
         'intervals in DOMAIN.'])
end

% Sort out the exponents:
exps = [];
if ( ~isempty(pref.singPrefs.exponents) )
    exps = pref.singPrefs.exponents;
    nExps = numel(pref.singPrefs.exponents);
    
    if ( nExps == 1 )
        
        % If only one exponent is supplied, assume the exponent at other
        % breakpoints are exactly same.
        exps = exps*ones(1,2*numIntervals);
        
    elseif ( nExps == 2 )
        
        % If the exponents are only supplied at endpoints of the entire
        % domain, then pad zeros at the interior breakpoints.
        exps = [exps(1) zeros(1, 2*(numIntervals-1)) exps(2)];
        
    elseif ( nExps == numIntervals + 1 )
        
        % If only one exponent is supplied for each interior breakpoint,
        % then we assume that the singularity take the same order on each
        % side.
        exps = exps(ceil(1:0.5:nExps - 0.5));
        
    elseif( nExps ~= 2*numIntervals )
        
        % The number of exponents supplied by user makes no sense.
        error('CHEBFUN:constructor', 'Invalid length for vector of exponents.');
    end
end

% Sort out the singularity types:
type = [];
if ( ~isempty(pref.singPrefs.singType) )
    type = pref.singPrefs.singType;
    nType = numel(pref.singPrefs.singType);
    
    if ( nType ~= 2*numIntervals )
        
        % If the number of exponents supplied by user isn't equal to twice the
        % the number of the FUNs, throw an error message:
        error('CHEBFUN:constructor', ['the number of the exponents is ' ...
            'inappropriate.']);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%  SPLITTING OFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In 'OFF' mode, seek only one piece with length < maxLength.
if ( ~pref.enableBreakpointDetection )
    % Set maximum length (i.e., number of sample points for CHEBTECH):
    maxn = pref.maxTotalLength;
    % Initialise the FUN array:
    funs{numIntervals} = fun.constructor();
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
        
        % Extract the exponents for this interval:
        if ( ~isempty(exps) )
            pref.singPrefs.exponents = exps(2*k-1:2*k);
        end
        
        % Replace the information for the singularity type in the preference:
        if ( ~isempty(type) )
            pref.singPrefs.singType = type(2*k-1:2*k);
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

% Set the maximum length (i.e., number of sample points for CHEBTECH):
pref.maxTotalLength = pref.breakpointPrefs.splitMaxTotalLength;
pref.techPrefs.maxLength = pref.breakpointPrefs.splitMaxLength;
% We extrapolate when splitting so that we can construct functions like
% chebfun(@sign,[-1 1]), which otherwise would not be happy at x = 0.
pref.techPrefs.extrapolate = true;

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
    
    % Extract the exponents for this interval:
    if ( ~isempty(exps) )
        pref.singPrefs.exponents = exps(2*k-1:2*k);
    end
    
    % Replace the information for the singularity type in the preference:
    if ( ~isempty(type) )
        pref.singPrefs.singType = type(2*k-1:2*k);
    end

    [funs{k}, ishappy(k), vscale] = ...
        getFun(opk, ends(k:k+1), vscale, hscale, pref);
    
    % For the case where vscale is Inf due to blowup in the interior of the
    % domain:
    if ( isinf(vscale) )
        % An infinite vscale doesn't help us at all, but ruin the consequential
        % constructions at the lower levels:
        vscale = 0; 
    end
    
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
    
    % Locate an edge/split location, compensating for exponents if necessary:
    if ( ~isempty(exps) && any( exps(2*k-1:2*k) ) )
        opkDetectEdge = @(x) opk(x)./((x - a).^exps(2*k - 1) .* ...
            (b - x).^exps(2*k));
        edge = chebfun.detectEdge(opkDetectEdge, [a, b], vscale, hscale);
    else
        edge = chebfun.detectEdge(opk, [a, b], vscale, hscale);
    end

    if ( ~isempty(exps) )
        % Before constructing the left FUN, sort out the exponents:
        pref.singPrefs.exponents = [exps(2*k-1) 0];
    end
    
    if ( ~isempty(type) )
        % Before constructing the left FUN, sort out the singType:
        pref.singPrefs.singType = [type(2*k-1) type(2*k-1)];
    end
    
    % Try to obtain happy child FUN objects on each new subinterval:
    [childLeft, happyLeft, vscale] = ...
        getFun(opk, [a, edge], vscale, hscale, pref);
    
    % For the case where vscale is Inf due to blowup in the interior of the
    % domain:
    if ( isinf(vscale) )
        % An infinite vscale doesn't help us at all, but ruin the consequential
        % constructions at the lower levels:
        vscale = 0;
    end
    
    if ( ~isempty(exps) )
        % Before constructing the right FUN, sort out the exponents:
        pref.singPrefs.exponents = [0 exps(2*k)];
    end
    
    if ( ~isempty(type) )
        % Before constructing the right FUN, sort out the singType:
        pref.singPrefs.singType = [type(2*k) type(2*k)];
    end
    
    [childRight, happyRight, vscale] = ...
        getFun(opk, [edge, b], vscale, hscale, pref);
    
    % For the case where vscale is Inf due to blowup in the interior of the
    % domain:
    if ( isinf(vscale) )
        % An infinite vscale doesn't help us at all, but ruin the consequential
        % constructions at the lower levels:
        vscale = 0;
    end
    
    % Check for happiness/sadness:
    sad = [sad(1:k-1), ~happyLeft, ~happyRight, sad(k+1:end)];

    % Insert new pieces in to existing funs:
    funs = [funs(1:k-1), {childLeft, childRight}, funs(k+1:end)];
    ends = [ends(1:k), edge, ends(k+1:end)];
    
    if ( ~isempty(exps) )
        exps = [exps(1:2*k-1), zeros(1,2), exps(2*k:end)];
    end
    
    if ( ~isempty(type) )
        type = [type(1:2*k-1), type(2*k-1), type(2*k) type(2*k:end)];
    end
    
    % If a cell was given, we must store pieces on new intervals:
    if ( iscell(op) )
        op = [op(1:k), {opk}, op(k+1:end)];
    end

    % Fail if too many points are required:
    if ( sum(cellfun(@length, funs)) > pref.maxTotalLength )
        warning('Function not resolved using %d pts.', ...
            sum(cellfun(@length, funs)));
        return
    end

end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GETFUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [g, ishappy, vscale] = getFun(op, domain, vscale, hscale, pref)
%GETFUN    Call the FUN constructor.

% If the interval is very small then skip adaptation and treat OP as a constant:
if ( diff(domain) < 4*1e-14*hscale && ~isnumeric(op) )
    op = op(mean(domain));
end

% Call the FUN constructor:
g = fun.constructor(op, domain, vscale, hscale, pref);
% See if the construction was happy:
ishappy = get(g, 'ishappy');
% Update the vertical scale:
vscale = max([vscale, get(g, 'vscale')]);

end
