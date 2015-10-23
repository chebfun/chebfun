function [f, mergedPts] = merge(f, index, pref)
%MERGE   Remove unnecessary breakpoints in from a CHEBFUN.
%   F = MERGE(F, PREF) removes unnecessary breakpoints from a CHEBFUN F. In
%   particular the kth breakpoint is removed if the resulting FUN on the
%   interval [x_{k-1}, x_{k+1}] can be represented to the same accuracy as F
%   with a fewer than PREF.MAXLENGTH points when PREF.SPLITTING = 0 and
%   PREF.SPLITPREFS.SPLITLENGTH points when PREF.SPLITTING = 1. If a PREF is
%   not passed, then the default CHEBFUN.PREF() is used.
%
%   [F, MERGEDPTS] = MERGE(F) returns the index of the merged endpoints in the
%   vector MERGEDPTS.
%
%   MERGE(F, INDEX) or MERGE(F, INDEX, PREF) attempts to eliminate the endpoints
%   specified in INDEX. MERGE(F, 'all') is equivalent to MERGE(F,
%   [2:length(F.domain)-1]). (Note that it doesn't make sense to consider
%   merging the first and final breakpoints.)
%
%   In all cases, elimination is attempted from left to right, and non-trivial
%   pointValues will prevent merging at the corresponding breakpoints.
%
%   Example:
%       f = chebfun(@(x) abs(x), 'splitting','on');
%       [g, mergedPts] = merge(f.^2);
%
% See also CHEBFUNPREF.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( numel(f) > 1 )
    % TODO:  Implement this.
    error('CHEBFUN:CHEBFUN:merge:quasi', ...
        'MERGE does not support quasimatrices.');
end

% Convert to a column CHEBFUN so that feval(f, x) returns a column instead of a
% row.  (Needed to construct a CHEBFUN out of @(x) feval(f, x) later on.)
isTrans = f.isTransposed;
if ( isTrans )
    f = f.';
end

% Parse the inputs:
if ( nargin == 1 )
    % Choose all indices by default:
    index = 2:numel(f.funs);
    % Obtain preferences:
    pref = chebfunpref();
elseif ( nargin == 2 )
    if ( isstruct(index) || isa(index, 'chebfunpref') ) % MERGE(F, PREF)
        % index actually is a struct of preferences or a CHEBFUNPREF
        pref = chebfunpref(index);
        % Choose all indices by default:
        index = 2:numel(f.funs);
    else                   % MERGE(F, INDEX)
        % indices passed, obtain preferences:
        pref = chebfunpref();
    end
end

[f, mergedPts] = mergeColumn(f, index, pref);

% Convert back to a row CHEBFUN if we started with one.
if ( isTrans )
    f.isTransposed = true;
end

end

function [f, mergedPts] = mergeColumn(f, index, pref)

% Deal with input arguments:
if ( isempty(index) )
    % No indices requested.
    mergedPts = [];
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
        % All the requested indices were trivial.
        mergedPts = [];
        return
    end

end

% Determine the maximum length of the merged pieces:
if ( ~pref.splitting )
    maxn = pref.techPrefs.maxLength;
else
    maxn = pref.splitPrefs.splitLength;
end
pref.techPrefs.maxLength = maxn;

% Splitting forces extrapolate:
if ( pref.splitting )
    pref.techPrefs.extrapolate = true;
end

% Obtain scales of the CHEBFUN:
vs = vscale(f);
hs = hscale(f);
pref.eps = max(eps, pref.eps);
mergedPts = [];

% Store data from input CHEBFUN:
oldPointVals = f.pointValues;
oldDom = f.domain;
oldFuns = f.funs;
newPointVals = oldPointVals;
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
    v = [oldPointVals(k,:); get(oldFuns{k-1}, 'rval'); get(oldFuns{k}, 'lval')];
    if ( all( norm(v([1, 1],:) - v(2:3,:), inf) >= 1e3*pref.eps ) || any(isinf(v(:))) )
        % Skip to next breakpoint:
        continue
    end

    % Call merge at the FUN level:
    [mergedFun, ishappy] = merge(newFuns{j-1}, newFuns{j}, vs, hs, pref);
    
    % Prevent merge if the result is not happy:
    if ( ~ishappy )
        % Skip to next breakpoint:
        continue
    end

    % Merging successful!
    mergedPts = [mergedPts, k]; %#ok<AGROW> % Store index of the removed point.
    newFuns{j-1} = mergedFun;   % Insert new fun.
    newFuns(j) = [];            % Remove unneeded fun.
    newDom = [newDom(1:j-1), newDom(j+1:end)]; % Update domain.
    % Update pointValues.
    newPointVals = [newPointVals(1:j-1,:) ; newPointVals(j+1:end,:)]; 

end

% Assign data to CHEBFUN:
f.domain = newDom;
f.funs = newFuns;
f.pointValues = newPointVals;

end

