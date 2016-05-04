function [funs, ends] = constructor(op, dom, data, pref)
%CONSTRUCTOR   CHEBFUN constructor.
%   FUNS = CONSTRUCTOR(OP, DOM) constructs the piecewise components (known as
%   "FUNS") used by a CHEBFUN object to represent the function OP on the
%   interval DOM. OP must be a function_handle, string, numerical vector, or a
%   cell array containing a combination of these first three data types. In the
%   later case, the number of elements in the array must be one less than the
%   length of the DOM vector.
%
%   It is not expected that CHEBFUN.CONSTRUCTOR() be called directly, 
%
%   If OP is a function_handle or a string, it should be vectorised in that it
%   accepts a column vector of length N and return a matrix of size N x M. If M
%   ~= 1, we say the resulting CHEBFUN is "array-valued".
%
%   CONSTRUCTOR(OP, DOM, DATA, PREF), where DATA is a MATLAB structure and PREF
%   is a CHEBFUNPREF object, allows construction data and alternative
%   construction preferences to be passed to the constructor.  See CHEBFUNPREF
%   for more details on preferences.
%
%   In particular, if PREF.SPLITTING = TRUE and OP is a function_handle or a
%   string, then the constructor adaptively introduces additional breakpoints
%   into the domain so as to better represent the function. These are returned
%   as the second output argument in [FUNS, END] = CONSTRUCTOR(OP, DOM).
%
%   The DATA structure input contains information which needs to be passed to
%   the lower layers about parameters which may affect the construction process.
%   Presently, the only fields CONSTRUCTOR expects DATA to have on input are
%   DATA.EXPONENTS and DATA.SINGTYPE, which convey information about endpoint
%   singularities. These fields are populated by CHEBFUN.PARSEINPUTS as need be.
%   Before calling the FUN constructor, DATA will be augmented to include
%   information about the construction domain as well as the horizontal and
%   vertical scales involved in the construction procedure.
%
% See also CHEBFUN, CHEBFUNPREF.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Initial setup:
numIntervals = numel(dom) - 1;

% Initialise hscale and vscale:
data.hscale = norm(dom, inf);
if ( isinf(data.hscale) )
    data.hscale = 1;
end
if ( isempty(data.vscale) )
    data.vscale = 0;
end

% Sanity check:
if ( iscell(op) && (numel(op) ~= numIntervals) )
    error('CHEBFUN:CHEBFUN:constructor:cellInput', ...
        ['Number of cell elements in OP must match the number of ', ...
         'intervals in DOMAIN.'])
end

% Construct the FUNs.
if ( pref.splitting )
    [funs, ends] = constructorSplit(op, dom, data, pref);
else
    [funs, ends] = constructorNoSplit(op, dom, data, pref);
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%  SPLITTING OFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [funs, ends] = ...
    constructorNoSplit(op, domain, data, pref)
% In 'OFF' mode, seek only one piece with length < maxLength.

% Initial setup:
numIntervals = numel(domain) - 1;
ends = domain;

% Initialise the FUN array:
funs = cell(1, numIntervals);

% We only want to throw the warning 'CHEBFUN:constructor:notResolved'once:
warningThrown = false;

singDetect = pref.blowup;
exps = data.exponents;
singTypes = data.singType;

% Loop over the intervals:
for k = 1:numIntervals
    endsk = ends(k:k+1);
    % Unwrap if OP is a cell:
    if ( iscell(op) )
        opk = op{k};
    else
        opk = op;
    end

    if ( singDetect )
        % Extract the exponents for this interval:
        if ( ~isempty(exps) )
            data.exponents = exps(2*k-1:2*k);
        end
        % Replace the information for the singularity type in the preference:
        if ( ~isempty(singTypes) && ~ischar(singTypes) )
            data.singType = singTypes(2*k-1:2*k);
        end
    end

    % Call GETFUN() (which calls FUN.CONSTRUCTOR()):
    [funs{k}, ishappy, data.vscale] = getFun(opk, endsk, data, pref);

    % Warn if unhappy (as we're unable to split the domain to improve):
    if ( ~ishappy && ~warningThrown )
        if ( strcmpi(func2str(pref.tech), 'trigtech') )
            str = 'a non-trig representation';
        else
            str = '''splitting on''';
        end
            
        warning('CHEBFUN:CHEBFUN:constructor:notResolved', ...
            ['Function not resolved using %d pts.', ...
            ' Have you tried ' str, '?'], pref.techPrefs.maxLength);
        warningThrown = true;
    end
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%  SPLITTING ON %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [funs, ends] = constructorSplit(op, dom, data, pref)
% In 'ON' mode, seek only many pieces with total length < maxlength.

% Initial setup:
numIntervals = numel(dom) - 1;
ends = dom;

% Set the maximum length (i.e., number of sample points for CHEBTECH):
pref.techPrefs.maxLength = pref.splitPrefs.splitLength;

% We extrapolate when splitting so that we can construct functions like
% chebfun(@sign,[-1 1]), which otherwise would not be happy at x = 0.
pref.techPrefs.extrapolate = true;

% Initialise the FUN array:
funs = cell(1, numIntervals);
% Initialise happiness:
ishappy = ones(1, numel(ends) - 1);

singDetect = pref.blowup;
exps = data.exponents;
singTypes = data.singType;

% Try to get one smooth piece for the entire domain before splitting:
for k = 1:numIntervals
    % Unwrap if OP is a cell:
    if ( iscell(op) )
        opk = op{k};
    else
        opk = op;
    end
    
    if ( singDetect )
        % Extract the singularity information for this interval:
        data = getSingInfo(exps, singTypes, 2*k-1:2*k, data);
    end

    [funs{k}, ishappy(k), data.vscale] = ...
        getFun(opk, ends(k:k+1), data, pref);
    
    % For the case where vscale is Inf due to blowup in the interior of the
    % domain:
    if ( isinf(data.vscale) )
        % An infinite vscale doesn't help us at all, but ruin the consequential
        % constructions at the lower levels:
        data.vscale = 0; 
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

    % Look for an edge:
    edge = fun.detectEdge(funs{k}, opk, data.hscale, data.vscale, pref);

    if ( singDetect )
        % Update singularity info:
        if ( ~isempty(exps) )
            exps = [exps(1:2*k-1), zeros(1,2), exps(2*k:end)];
        end
        if ( ~isempty(singTypes) )
            singTypes = [singTypes(1:2*k-1), singTypes(2*k-1), ...
                singTypes(2*k) singTypes(2*k:end)];
        end
    end
    
    if ( singDetect )
        % Extract the singularity information for this interval:
        data = getSingInfo(exps, singTypes, 2*k-1:2*k, data);
    end
    
    % Try to obtain happy child FUN objects on each new subinterval:
    [childLeft, happyLeft, data.vscale] = getFun(opk, [a, edge], data, pref);
    
    if ( singDetect )
        % Extract the singularity information for this interval:
        data = getSingInfo(exps, singTypes, 2*k+1:2*k+2, data);
    end
    
    [childRight, happyRight, data.vscale] = getFun(opk, [edge, b], data, pref);
    
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
    len = sum(cellfun(@length, funs));
    if ( len > pref.splitPrefs.splitMaxLength )
        warning('CHEBFUN:CHEBFUN:constructor:funNotResolved', ...
            'Function not resolved using %d pts.', ...
            sum(cellfun(@length, funs)));
        return
    end

end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GETFUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [g, ishappy, vscale] = getFun(op, dom, data, pref)
%GETFUN   Call the FUN constructor.

% If the interval is very small then skip adaptation and treat OP as a constant:
if ( diff(dom) < 4*1e-14*data.hscale && ~isnumeric(op) )
    op = op(mean(dom));
end

% Bolt domain information onto the data structure.
data.domain = dom;

% Call the FUN constructor:
g = fun.constructor(op, data, pref);

% See if the construction was happy:
ishappy = get(g, 'ishappy');

% Update the vertical scale:
if ( ishappy )
    vscale = max([data.vscale, get(g, 'vscale')]);
else
    vscale = data.vscale;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GETSINGINFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = getSingInfo(exps, singTypes, kk, data)
% Place information about the singularity type in the data structure.
if ( ~isempty(exps) )
    data.exponents = exps(kk);
end
if ( ~isempty(singTypes) )
    data.singType = singTypes(kk);
end
end
