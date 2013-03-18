function f = restrict(f, s)
%RESTRICT Restrict a CHEBTECH to a subinterval.
%   RESCTRICT(F, S) returns a CHEBTECH that is restricted to the subinterval
%   [S(1),S(2)] of [-1, 1]. Note that that since CHEBTECH only live on [-1,1], a
%   linear change of variables is implicitly applied.
%
%   If length(S) > 2, i.e., S = [S1, S2, S3, ...], then RESCTRICT(F, S) returns
%   an array of CHEBTECH objects, where the entries hold F restricted to each of
%   the subintervals defined by S.
%
%   If F is a vector-valued function, say [F1, F2], then the restrict(F, S =
%   [S1, S2, S3]) returns the vector-valued CHEBTECH {restrict(F1,S).
%   restrict(F2, S)}.
%
%   Note that restrict does not 'simplify' its output.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    return
end

% Check if subint is actually a subinterval
if ( s(1) < -1 || s(end) > 1 || any(diff(s) <= 0) )
    error('CHEBTECH:restrict:badinterval', 'Not a valid interval.')
elseif ( numel(s) == 2 && all(s == [-1, 1]) )
    % Nothing to do here!
    return
end

% Grab sizes:
[n, m] = size(f);
numInts = numel(s) - 1;

% Compute new values on the grid:
x = f.chebpts(n);                            % old grid
y = .5*[1-x, 1+x] * [s(1:end-1) ; s(2:end)]; % new grid
values = feval(f, y);                        % new values

% If F is multivalued, we must rearrange the order of the columns:
% (e.g., [a1 a2 b1 b2 c1 c2] -> [a1 b1 c1 a2 b2 c3] => idx = [1 3 5 2 4 6].
if ( m > 1 )
    numCols = m*numInts;
    idx = reshape(reshape(1:numCols, numInts, m)', 1, numCols);
    values = values(:, idx);
end

% Update coeffs and vscale:
coeffs = f.chebpoly(values);                 % new coeffs
vscale = max(abs(values), [], 1);            % new vscale

% Append data to fun:
f.values = values;
f.coeffs = coeffs;
f.vscale = vscale;
% f.epslevel = f.epslevel; % epslevel does not change (or become a vector yet).

% Convert to an array:
f = mat2cell(f, repmat(m, 1, numInts));

end

