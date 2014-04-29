function [f, mergedPts] = merge(f, index, pref)
%MERGE   Remove unnecessary breakpoints in from a CHEBFUN.
%   F = MERGE(F, PREF) removes unnecessary breakpoints from a CHEBFUN F. In
%   particular the kth breakpoint is removed if the resulting FUN on the
%   interval [x_{k-1}, x_{k+1}] can be represented with a fewer than
%   PREF.MAXLENGTH points when PREF.ENABLEBREAKPOINTDETECTION = 0 and
%   PREF.BREAKPOINTPREFS.SPLITMAXLENGTH points when
%   PREF.ENABLEBREAKPOINTDETECTION = 1. If a PREF is not passed, then the
%   default CHEBFUN.PREF() is used.
%
%   [F, MERGEDPTS] = MERGE(F) returns the index of the merged endpoints in the
%   vector MERGEDPTS.
%
%   MERGE(F, INDEX) attempts to eliminate the endpoints specified in INDEX.
%   MERGE(F, 'all') is equivalent to MERGE(F, [2:length(F.domain)-1]). (Note
%   that it doesn't make sense to consider merging the first and final
%   breakpoints.)
%
%   In all cases, elimination is attempted from left to right, and non-trivial
%   pointValues will prevent merging at the corresponding breakpoints.
%
%   Example:
%       f = chebfun(@(x) abs(x), 'splitting','on');
%       [g, mergedPts] = merge(f.^2);
%
% See also CHEBPREF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

if ( numel(f) > 1 )
    error('CHEBFUN:merge:quasi', 'MERGE does not support quasimatrices.');
end

% Parse the inputs:
if ( nargin == 1 )
    % Choose all indices by default:
    index = 2:numel(f.funs);
    % Obtain preferences:
    pref = chebpref();
elseif ( nargin == 2 )
    if ( isstruct(index) || isa(index, 'chebpref') ) % MERGE(F, PREF)
        % index actually is a struct of perefernces or a CHEBPREF
        pref = chebpref(index);
        % Choose all indices by default:
        index = 2:numel(f.funs);
    else                   % MERGE(F, INDEX)
        % indices passed, obtain preferences:
        pref = chebpref();
    end
end

% Deal with input arguments:
if ( isempty(index) )
    % No indices requested!
    return

elseif ( ischar(index) )
    % 'all' indices requested:
    index = 2:numel(f.funs);

else
    % Index of endpoints was provided:
    index = unique(index);
    % Break points must be in range 2:numel(f.funs)
    index((index <= 1) | (index > numel(f.funs))) = [];
    
    if ( isempty(index) )
        % All the requested indices were trivial:
        return
    end

end

% Determine the maximum length of the merged pieces:
if ( ~pref.enableBreakpointDetection )
    maxn = pref.techPrefs.maxLength;
else
    maxn = pref.breakpointPrefs.splitMaxLength;
end
pref.techPrefs.maxLength = maxn;

% Splitting forces extrapolate:
if ( pref.enableBreakpointDetection )
    pref.techPrefs.extrapolate = true;
end

% Obtain scales of the CHEBFUN:
vs = vscale(f);
hs = hscale(f);
tol = epslevel(f);
mergedPts = [];

% Store data from input CHEBFUN:
oldPointValues = f.pointValues;
oldDom = f.domain;
oldFuns = f.funs;
newPointValues = oldPointValues;
newDom = oldDom;
newFuns = oldFuns;

% Loop through the index:
for k = index
    % Find corresponding break:
    j = find(oldDom(k) == newDom, 1, 'first');
    % And lengths of funs on either side:
    lengthPrevFun = length(newFuns{j-1});
    lengthThisFun = length(newFuns{j});

    % Prevent merge if existing FUN lengths add to more than 1.2*maxn:
    if ( lengthPrevFun + lengthThisFun >= 1.2*maxn )
        % Skip to next breakpoint:
        continue
    end

    % Prevent merging if there are jumps:
    v = [ oldPointValues(k,:) ; get(oldFuns{k-1}, 'rval') ; get(oldFuns{k}, 'lval') ];
    if ( norm(v([1,1],:) - v(2:3,:), inf) >= 1e3*tol )
        % Skip to next breakpoint:
        continue
    end

    % Create copy of f with values at endpoints of interval on which we are
    % trying to create a merged fun replaced by left- and right-hand limits.
    % If the underlying tech samples right at the endpoints, this prevents
    % discontinuities that occur there from messing things up.
    g = f;
    g.pointValues(k-1,:) = get(oldFuns{k-1}, 'lval');
    g.pointValues(k+1,:) = get(oldFuns{k}, 'rval');

    % Grab the correct exponents:
    
    if ( issing(f) )
        if ( ~issing(newFuns{j-1}) && ~issing(newFuns{j}) )
            pref.singPrefs.exponents = [];
        elseif ( issing(newFuns{j-1}) && ~issing(newFuns{j}) )
            exps = get(newFuns{j-1}, 'exponents');
            pref.singPrefs.exponents = [exps(1) 0];
        elseif ( ~issing(newFuns{j-1}) && issing(newFuns{j}) )
            exps = get(newFuns{j}, 'exponents');
            pref.singPrefs.exponents = [0 exps(2)];
        else
            expsLeft = get(newFuns{j-1}, 'exponents');
            expsRight = get(newFuns{j}, 'exponents');
            pref.singPrefs.exponents = [expsLeft(1) expsRight(2)];
        end
    end
    
    % Attempt to form a merged FUN:
    mergedFun = fun.constructor(@(x) feval(g, x),  ...
        [newDom(j-1), newDom(j+1)], vs, hs, pref);

    % Prevent merge if the result is not happy:
    if ( ~get(mergedFun, 'ishappy') )
        % Skip to next breakpoint:
        continue
    end

    % Merging successful!
    mergedPts = [mergedPts, k]; %#ok<AGROW> % Store index of the removed point.
    newFuns{j-1} = mergedFun;   % Insert new fun.
    newFuns(j) = [];            % Remove unneeded fun.
    newDom = [newDom(1:j-1), newDom(j+1:end)];        % Update domain.
    newPointValues = [newPointValues(1:j-1,:) ; newPointValues(j+1:end,:)]; % Update pointValues.

end

% Assign data to CHEBFUN:
f.domain = newDom;
f.funs = newFuns;
f.pointValues = newPointValues;

end
