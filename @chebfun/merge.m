function [f, mergedPts] = merge(f, index, pref)
%MERGE   Remove unnecessary breakpoints in from a CHEBFUN.
%   F = MERGE(F, PREF) removes unnecessary breakpoints from a CHEBFUN F. In
%   particular the kth breakpoint is removed if the resulting FUN on the
%   interval [x_{k-1}, x_{k+1}] can be represented with a fewer than
%   PREF.CHEBFUN.MAXDEGREE points when PREF.CHEBFUN.SPLITTING = 0 and
%   PREF.CHEBFUN.SPLITDEGREE points when PREF.CHEBFUN.SPLITTING = 1. If a PREF
%   is not passed, then the default CHEBFUN.PREF() is used.
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
%   impulses will prevent merging at the corresponding breakpoints.
%
%   Example:
%       f = chebfun(@(x) abs(x), 'splitting','on');
%       [g, mergedPts] = merge(f.^2);
%
% See also SPLITTING, PREF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Parse the inputs:
if ( nargin == 1 )
    % Choose all indices by default:
    index = 2:numel(f.funs);
    % Obtain preferences:
    pref = chebfun.pref();
elseif ( nargin == 2 )
    if ( isstruct(index) ) % MERGE(F, PREF)
        pref = index;
        % Choose all indices by default:
        index = 2:numel(f.funs);
    else                   % MERGE(F, INDEX)
        % Obtain preferences:
        pref = chebfun.pref();
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
    % Break points must be in in range 2:numel(f.funs)
    index((index <= 1) | (index > numel(f.funs))) = [];

end

% Determine the maximum length of the merged pieces:
if ( ~pref.chebfun.splitting )
    maxn = pref.chebfun.maxdegree + 1;
else
    maxn = pref.chebfun.splitdegree + 1;
end
pref.chebfun.maxSamples = maxn;

% Obtain scales of the CHEBFUN:
vs = vscale(f);
hs = hscale(f);
tol = epslevel(f);
mergedPts = [];

% Store data from input CHEBFUN:
oldImps = f.impulses;
oldDom = f.domain;
oldFuns = f.funs;
newImps = oldImps;
newDom = oldDom;
newFuns = oldFuns;

% Merge preferences for FUN constructor call:
pref = fun.pref(pref, pref.chebfun);

% Loop through the index:
for k = index

    % Find corresponding break:
    j = find(oldDom(k) == newDom, 1, 'first');
    % And lengths of funs on either side:
    ljm1 = length(newFuns{j-1});
    lj1 = length(newFuns{j});

    % Prevent merge if existing FUN lengths add to more than 1.2*maxn:
    if ( ljm1 + lj1 >= 1.2*maxn )
        % Skip to next breakpoint:
        continue
    end

    % Prevent merge if nontrivial impulses:
    if ( any(any(oldImps(k, :, 2:end))) )
        % Skip to next breakpoint:
        continue
    end

    % Prevent merging if there are jumps:
    v = [ oldImps(k,:,1), get(oldFuns{k-1}, 'rval'), get(oldFuns{k}, 'lval') ];
    if ( norm(v(1) - v(2:3), inf) >= 1e3*tol )
        % Skip to next breakpoint:
        continue
    end

    % Attempt to form a merged FUN:
    mergedFun = fun.constructor(@(x) feval(f, x),  ...
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
    newImps = [newImps(1:j-1,:,:) ; newImps(j+1:end,:,:)]; % Update impulses.

end

% Assign data to CHEBFUN:
f.domain = newDom;
f.funs = newFuns;
f.impulses = newImps;

end
