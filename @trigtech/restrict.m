function F = restrict(f, s)
%RESTRICT   Restrict a TRIGTECH to a subinterval.
%   RESTRICT(F, S) returns a TRIGTECH that is restricted to the subinterval
%   [S(1),S(2)] of [-1,1]. Note that that since TRIGTECH objects only live on
%   [-1,1], a linear change of variables is implicitly applied.
%
%   If length(S) > 2, i.e., S = [S1, S2, S3, ...], then RESTRICT(F, S) returns
%   a cell-array of TRIGTECH objects, where the entries hold F restricted to
%   each of the subintervals defined by S.
%
%   If F is an array-valued function, say [F1, F2], then the restrict(F, S =
%   [S1, S2, S3]) returns the array-valued TRIGTECH {restrict(F1,S).
%   restrict(F2, S)}.
%
%   Note that restrict does not 'simplify' its output.
%
%   Warning: If F is not also smooth and periodic on S, then the resulting
%   TRIGTECH will not be happy.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    F = f;
    return
end

% Check if s is actually a subinterval:
if ( (s(1) < -1) || (s(end) > 1) || (any(diff(s) <= 0)) )
    error('TRIGTECH:restrict:badinterval', 'Not a valid interval.')
elseif ( (numel(s) == 2) && all(s == [-1, 1]) )
    % Nothing to do here!
    F = f;
    return
end

% Determine number of intervals.
numInts = numel(s) - 1;

if ( numInts > 1 )
    error('CHEBFUN:TRIGTECH:restrict:multIntervals', ...
        'Cannot restrict a TRIGTECH to multiple intervals.');
end

op = @(x) feval(f, .5 * [1 - x, 1 + x] * [s(1) ; s(end)]);
data.vscale = vscale(f);
pref = f.techPref;
pref.minSamples = min(length(f), pref.minSamples);
F = f.make(op, data, pref);

if ( f.ishappy && ~F.ishappy )
    error('CHEBFUN:TRIGTECH:restrict:notPeriodic', ...
        'Restrict failed. Perhaps f is not periodic in [s(1), s(2)]?');
end

end
