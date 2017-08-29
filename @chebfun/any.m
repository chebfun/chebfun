function a = any(f, dim)
%ANY   True if any value of a CHEBFUN is nonzero. ANY ignores entries that are
%      NaN (Not a Number).
%   ANY(X, DIM), where X is an array-valued CHEBFUN, works down the dimension
%   DIM. If DIM is the CHEBFUN (continuous) dimension, then ANY returns a
%   logical column vector (or row) in which the Jth element is TRUE if the Jth
%   column (or row) has a nonzero value. Otherwise, ANY returns a CHEBFUN which
%   takes the value 1 wherever any of the columns (or rows) of X are nonzero,
%   and zero everywhere else.
%
%   ANY(X) is shorthand for ANY(X, 1).
%
% See also ALL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information

% Parse inputs:
if ( nargin < 2 )
    dim = 1;
end

if ( (dim ~= 1) && (dim ~= 2) )
    error('CHEBFUN:CHEBFUN:any:dim', 'DIM input must be 1 or 2.');
end

% Deal with row CHEBFUNs by transposing first:
isTransposed = f.isTransposed;
if ( isTransposed )
    f = f.';
    dim = mod(dim, 2) + 1; % Maps 1 to 2 and 2 to 1.
end

if ( dim == 1 )
    a = anyDim1(f);     % ANY() down the columns (continuous dimension).
else
    a = anyDim2(f);     % ANY() across the rows (discrete dimension).
end

% Transpose output, if necessary:
if ( isTransposed )
    a = a .';
end

end

function a = anyDim1(f)
%ANYDIM1   Column-wise ANY() for a column CHEBFUN.

% NB:  This code handles empty f without an extra check, since in this case,
% f.pointValues is empty and numel(f.funs) = 0, so the loop will never execute.

% Check the pointValues:
a = any(f.pointValues);

% If necessary, check all of the FUNs:
k = 1;
while ( ~all(a) && (k <= numel(f.funs)) )
    a = a | any(f.funs{k});
    k = k + 1;
end

end

function a = anyDim2(f)
%ANYDIM2   Row-wise ANY() for a column CHEBFUN.

% Deal with the empty case:
if ( isempty(f) )
    a = f;
    return
end

% Zero small breakpoint values so they will be seen as zero by ANY():
f = thresholdBreakpointValues(f);

% Split up f at its roots and call ANY() for the individual FUNs:
a = addBreaksAtRoots(f);
for k = 1:1:numel(a.funs)
    a.funs{k} = any(a.funs{k}, 2);
end

% Deal with the pointValues:
% [TODO]; not sure about this:
% Deal with the impulses:
%   1.  Call ANY() along the third dimension to pick up higher-order impulses.
% a.impulses = any(any(a.impulses, 3), 2);

a.pointValues = any(any(a.pointValues, 3), 2);

% Get rid of unnecessary breakpoints:
a = merge(a);

end
