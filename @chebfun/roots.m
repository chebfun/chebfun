function r = roots(F, varargin)
%ROOTS   Roots of a CHEBFUN.
%   ROOTS(F) returns the roots of F in its domain of definition. By default,
%   roots are returned at jumps in F which pass through zero, and if F is
%   identically zero on a part of its domain, then a single root is returned at
%   the midpoint. Each of these behaviours can be modified using the optional
%   inputs described below:
%
%   ROOTS(F, 'nojump') prevents ROOTS() from returning points where F changes
%   sign due to a jump discontinuity, such as roots(chebfun(@sign, 'splitting',
%   'on')).
%
%   ROOTS(F, 'nozerofun') prevents ROOTS() from returning a zero at the midpoint
%   of the domain of F when F if identically zero, such as ROOTS(chebfun(0)).
% 
%   ROOTS(F, 'norecursion') deactivates the recursion procedure used to compute
%   roots (see the Guide 3: Rootfinding and minima and maxima for more
%   information on this recursion procedure).
%
%   ROOTS(F, 'all') returns the roots of all the FUN objects representing the
%   smooth pieces of F. Note that by default this disables recursion, and so is
%   equivalent to ROOTS(F, 'all', 'norecursion').
%
%   ROOTS(F, 'complex') returns the real and complex roots of the FUN objects
%   representing the smooth pieces of F that are determined to be non-spurious.
%   (See CHEBELLIPSEPLOT). This capability may remove some spurious roots that
%   can appear if using ROOTS(F, 'all'). ROOTS(F, 'complex') is equivalent to
%   ROOTS(F, 'complex', 'recursion').
%
%   ROOTS(F, 'all', 'recursion') and ROOTS(F, 'complex','norecursion') can be
%   used to activate and deactivate the recursion procedure respectively, to
%   compute the roots as explained in the 'all' and 'complex' modes.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: Scales and tolerances are quite arbitrary here..

% TODO: Document array-valued and quasimatrix cases.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note: 
%  Here we get around the fact that FUN/ROOTS will return NaN values for
%  array-valued FUN objects by simply ignoring NaNs whenever we come across
%  them. In particular, we make use of the fact that the built in MAX() command
%  ignores NaNs. Once we have looped through all the intervals, we sort the
%  resulting matrix of roots (note that SORT() also treats NaN as greater than
%  any double, whereas MAX() in effect treats it as less than any double) and
%  then remove any rows which contain only NaN values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Deal with the trivial empty case:
if ( isempty(F) )
    return
end

rootsPref = parseInputs(F, varargin{:});

% Loops over the columns:
r = cell(1, numel(F));
for k = 1:numel(F)
    r{k} = columnRoots(F(k), rootsPref); 
end

% Pad with NaNs for quasimatrices:
if ( numel(F) > 1 )
    l = cellfun(@length, r);
    ml = max(l);
    for k = 1:numel(F)
        r{k} = [r{k} ; NaN(ml-l(k), 1)];
    end
end

r = cell2mat(r);

end

function r = columnRoots(f, rootsPref)

% Set horizontal and vertical scales:
el = epslevel(f);
hs = hscale(f);
vs = vscale(f);
htol = 100*eps*hs;
vtol = el*vs;
dom = f.domain;

% Initialise vector to store roots:
numCols = size(f.funs{1}, 2);
NaNRow = NaN(1, numCols);
r = NaNRow;

% Zero pointValues are roots.
index = abs(f.pointValues(1,:)) < vs*el;
if ( any(index) )
    % Left pointValues is zero: (or sufficiently close)
    r(index) = dom(1);
end

funs = f.funs;
numFuns = numel(funs);
for k = 1:numFuns

    %% Roots within the subdomains:
    % Get the roots of the current fun:
    rk = roots(funs{k}, rootsPref);
    % Trim out roots that are repeated on either side of the breakpoint:
    if ( ~isempty(r) )
        % NOTE: Here we use the fact that MAX() ignores NaN values!
        rk(abs(bsxfun(@minus, max(r, [], 1), rk)) < htol) = NaN;
    end
    % Append new roots to r:
    r = [ r ; rk ]; %#ok<AGROW>
        
    %% Look for roots at next breakpoint:
    index = false(1, numCols);
    if ( rootsPref.breakRoot )
        index = abs(f.pointValues(k+1,:)) < vtol; % Include if zero pointValues
    end
    if ( rootsPref.jumpRoot && k < numFuns )     % Or a change in sign in a jump
        index = index | ( get(funs{k}, 'rval').*get(funs{k+1}, 'lval') <= 0 );
    end
    if ( ~isempty(r) )                           % But not if already a root!
        index = index & ~(abs(max(r, [], 1) - dom(k+1)) < htol);
    end
    rk = NaNRow;
    rk(index) = dom(k+1);
    % Append new roots to r:
    r = [ r ; rk ]; %#ok<AGROW>    

end

% Remove unnecessary NaNs:
r = sort(r, 1);             % Sort will place NaNs in the final rows.
r(all(isnan(r), 2),:) = []; % This removes any rows which contain only NaNs.

end

function rootsPref = parseInputs(f, varargin)
% Parse the preferences. See documentation in ROOTS().

% Defaults:
rootsPref = struct('all', 0, 'recurse', 1, 'prune', 0,  'zeroFun', 1, ...
    'jumpRoot', 1, 'breakRoot', 1, 'qz', 0, 'filter', []);

if ( ~isreal(f) )
    % 'jumpRoots' only makes sense for real-valued functions, so disable it:
    rootsPref.jumpRoot = false;
end

% Parse the inputs:
recurseHasBeenSet = 0;
for k = 1:numel(varargin)
    argin = lower(varargin{k});
    switch ( argin )
        case 'all'
            rootsPref.all = 1;
            rootsPref.prune = 0;
            if ( ~recurseHasBeenSet )
                rootsPref.recurse = 0;
            end
        case 'complex'
            rootsPref.prune = 1;
            rootsPref.all = 1;
        case 'zerofun'
            rootsPref.zeroFun = 1;
        case 'nozerofun'
            rootsPref.zeroFun = 0;
        case 'jump'
            rootsPref.jumpRoot = 1;
        case 'nojump'
            rootsPref.jumpRoot = 0;
        case 'breaks'
            rootsPref.breakRoot = 1;
        case 'nobreaks'
            rootsPref.breakRoot = 0;
        case 'qz'
            rootsPref.qz = 1; 
        case {'recursion', 'recurse'}
            rootsPref.recurse = 1;
            recurseHasBeenSet = 1;
        case {'norecursion', 'norecurse'}
            rootsPref.recurse = 0;
        otherwise
            error('CHEBFUN:CHEBFUN:roots:parseInputs:UnknownOption', ...
                'Unknown option in ROOTS.')
    end
end

end
