function [y, x] = minandmax(f, flag)
%MINANDMAX   Minimum and maximum values of a CHEBFUN.
%   Y = MINANDMAX(F) returns the range of the CHEBFUN F such that Y(1,1) =
%   min(F) and Y(2,1) = max(F).
%
%   [Y, X] = MINANDMAX(F) returns also points X such that F(X(j,1)) = Y(j,1), j
%   = 1, 2.
%
%   [Y, X] = MINANDMAX(F, 'local') returns not just the global minimum and
%   maximum values, but all of the local extrema (i.e., local min and max).
%   Note that point values are not regarded as local extrema.
%
%   If F is complex-valued, absolute values are taken to determine extrema, but
%   the resulting values correspond to those of the original function.
%   
%   If F is array-valued, then the columns of X and Y correspond to the columns
%   of F. NaNs are used to pad Y and X when the 'local' flag is used and the
%   columns are not of the same length.
%
% See also MAX, MIN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Trivial empty case
if ( isempty(f) )
    x = [];
    y = [];
    return
end

% Number of columns of an array-valued CHEBFUN:
numCols = size(f.funs{1}, 2);

if ( (nargin > 1) && strcmpi(flag, 'local') )
    % Deal with local case:
    [y, x] = localMinAndMax(f);
    return
end

% Deal with array-valued case:
if ( numCols > 1 )
    f = mat2cell(f);
    x = zeros(2, numCols);
    y = zeros(2, numCols);
    for k = 1:numel(f)
        [y(:,k), x(:,k)] = minandmax(f{k});
    end
    return
end

% Deal with complex case:
if ( ~isreal(f) )
    realf = real(f);
    imagf = imag(f);
    g = realf.*realf + imagf.*imagf;
    g = simplify(g);
    [ignored, x] = minandmax(g);
    y = feval(f, x);
    return
end

% NOTE: From here onwards, f will only be a scalar-valued CHEBFUN.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPULSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = [inf, -inf];
x = [NaN, NaN];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SMOOTH PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dom = f.domain;
nfuns = numel(f.funs);
yy = [zeros(nfuns, 2) ; y];
xx = [zeros(nfuns, 2) ; x];
for k = 1:nfuns
    [yy(k,:), xx(k,:)] = minandmax(f.funs{k});
end
[y(1), I1] = min(yy(:,1));
[y(2), I2] = max(yy(:,2));

x(1) = xx(I1,1);
x(2) = xx(I2,2);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%% BREAKPOINT VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind = find(f.pointValues < y(1));
if ( ~isempty(ind) )
    [y(1), k] = min(f.pointValues(ind,1));
    x(1) = dom(ind(k));
end
ind = find(f.pointValues(:,1) > y(2));
if ~isempty(ind)
    [y(2), k] = max(f.pointValues(ind,1));
    x(2) = dom(ind(k));
end

% Output column vector:
y = y.';
x = x.';

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCAL EXTREMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y, x] = localMinAndMax(f)
%LOCALMINANDMAX   Compute local extrema of f.

% Compute the turning points:
df = diff(f);
% Ensure endpoints are included:
df.pointValues([1,end],:) = 0;
% Call ROOTS():
x = roots(df);

% Evaluate the function at the turning points:
y = zeros(size(x));
for k = 1:size(x, 2)
    % Evaluate at the kth column in x:
    tmp = feval(f, x(:,k));
    % Extract the kth column in the output:
    y(:,k) = tmp(:,k);
end
% Ensure f(NaN) = NaN for array-valued CHEBFUN objects:
y(isnan(x)) = NaN;

end
