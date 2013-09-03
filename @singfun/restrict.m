function g = restrict(f, s)
%RESTRICT   Restrict a SINGFUN to a subinterval.
%   RESCTRICT(F, S) returns a SINGFUN that is restricted to the subinterval
%   [S(1), S(2)] of [-1,1]. Note that since SINGFUN objects only live on
%   [-1,1], a linear change of variables is implicitly applied.
%
%   If length(S) > 2, i.e., S = [S1, S2, S3, ...], then RESCTRICT(F, S) returns
%   an array of SINGFUN objects, where the entries hold F restricted to each of
%   the subintervals defined by S.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    return
end

% Check if s is actually a subinterval:
if ( (s(1) < -1) || (s(end) > 1) || (any(diff(s) <= 0)) )
    error('CHEBFUN:SINGFUN:restrict:badinterval', 'Not a valid interval.')
elseif ( (numel(s) == 2) && all(s == [-1, 1]) )
    % Nothing to do here!
    return
end

% Number of subintervals:
numInts = numel(s) - 1;

% Preallocate a cell
g = cell(1, numInts);

% Set the bound which is deemed 'finite':
fntbnd = realmax;

% Absolute values at the endpoints of the subintervals:
endVal = zeros(1, numInts+1);
for j = 1:numInts+1
    endVal(j) = abs(feval(f, s(j)));
end

% In the following construction process, we try to make the best of the
% information stored in the original singfun f.

for j = 1:numInts
    % Check if the first subinterval includes the left endpoint of the original
    % domain, i.e. -1 and if the function value at the right endpoint of the
    % last subinterval is finite. This rules out the possibility that the right
    % endpoint of the last subinterval is too close to 1, when a pole is
    % present at 1.
    if ( (s(j) == -1) && (endVal(j + 1) < fntbnd) )
        % define the new operator
        op = @(x) feval(f, ((1 - x)*s(1) + (1 + x)*s(2))/2);

        % call the singfun constructor
        gtmp = singfun( op, [f.exponents(1) 0], {'sing', 'none'}, [] );
        g{1} = gtmp;

        continue
    end

    % Check if the last subinterval includes the right endpoint of the original
    % domain, i.e. 1 and if the function value at the left endpoint of the
    % first subinterval is finite. This rules out the possibility that the left
    % endpoint of the first subinterval is too close to -1, when a pole is
    % present at -1.
    if ( (endVal(j) < fntbnd) && (s(j+1) == 1) )
        % define the new operator
        op = @(x) feval(f, ((1 - x)*s(end-1) + (1 + x)*s(end))/2);

        % call the singfun constructor
        gtmp = singfun( op, [0 f.exponents(2)], {'none', 'sing'}, [] );

        % put in cell
        g{end} = gtmp;

        continue
    end

    % Check if any of s(1) and s(end) is far from a singularity at -1 or 1.
    if ( max(endVal(j:j+1)) < fntbnd )
        % define the new operator which will be evaluated in the subinterval by
        % the singfun ctor.
        op = @(x) feval(f, ((1 - x)*s(j) + (1 + x)*s(j+1))/2);

        % call the singfun constructor
        gtmp = singfun( op, zeros(1, 2), {'none', 'none'}, [] );

        % put in cell
        g{j} = gtmp;

        continue
    end

    % For all remaining cases, call the constructor. This is mainly used for
    % handling subintervals with one of the endpoints too close to a
    % singularity.

    % define the new operator
    op = @(x) feval(f, ((1 - x)*s(j) + (1 + x)*s(j+1))/2);

    % call the singfun constructor
    gtmp = singfun( op, [], {'sing', 'sing'}, [] );

    % put in cell
    g{j} = gtmp;
end

% If there is only one subinterval, return the singfun, instead of a cell.
if ( numInts == 1 )
    g = g{:};
end

end
