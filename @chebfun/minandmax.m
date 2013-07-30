function [y, x] = minandmax(f, flag)
% MINANDMAX Minimum and maximum values of a chebfun.
%
% Y = MINANDMAX(F) returns the range of the chebfun F such that Y(1) = min(F)
% and Y(2) = max(F).
%
% [Y, X] = MINANDMAX(F) returns also points X such that F(X(j)) = Y(j), j = 1,2.
%
% [Y, X] = MINANDMAX(F, 'local') returns not just the global minimum and maximum
% values, but all of the local extrema (i.e. local min and max).
%
% If F is complex-valued, absolute values are taken to determine extrema, but
% the resulting values correspond to those of the original function.
%
% See also CHEBFUN/MAX, CHEBFUN/MIN.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with complex case:
if ( nargin > 1 && strcmpi(flag, 'local') )
    [y, x] = localminandmax(f);
    return
end

% [TODO]: Deal with complex case:
% if ( ~isreal(f) )
%     realf = real(f);
%     imagf = imag(f);
%     h = realf.*realf + imagf.*imagf;
%     h = simplify(h);
%     [y, x] = minandmax(h);
%     y = feval(f, x);
%     y = y(:,1:size(x,2)+1:end);
%     return
% end

y = [inf, -inf];
x = [inf, inf];

% Negative impulse, return y(1) = -inf
ind = find(min(f.impulses(2:end,:), [], 1) < 0 ,1 , 'first');
if ( ~isempty(ind) )
    y(1) = -inf;
    x(1) = f.domain(ind);
end

% Positive impulse, return y(2) = -inf
ind = find(max(f.impulses(2:end,:), [], 1) > 0, 1, 'first');
if ( ~isempty(ind) )
    y(2) = inf;
    x(2) = f.domain(ind);
end

if ( all(isfinite(x)) )
    % We're done.
    return
end

dom = f.domain;
nfuns = numel(f.funs);
yy = [zeros(nfuns,2) ; y];
xx = [zeros(nfuns,2) ; x];
for k = 1:nfuns
    [yy(k,:), xx(k,:)] = minandmax(f.funs{k});
end
[y(1), I1] = min(yy(:,1));
[y(2), I2] = max(yy(:,2));

x(1) = xx(I1,1);
x(2) = xx(I2,2);

% Check values at end points:
ind = find(f.impulses(1,:) < y(1));
if ( ~isempty(ind) )
    [y(1), k] = min(f.impulses(1,ind));
    x(1) = dom(ind(k));
end
ind = find(f.impulses(1,:) > y(2));
if ~isempty(ind)
    [y(2), k] = max(f.impulses(1,ind));
    x(2) = dom(ind(k));
end

end

function [y, x] = localminandmax(f)

% Compute the turning points:
df = diff(f);
x = roots(df);

% Deal with the end-points:
dom = f.domain;
xends = dom([1, end]);
dfends = feval(df, xends);
if ( dfends(1)~=0 )
    x = [xends(1) ; x];
end
if ( dfends(2)~=0 )
    x = [x ; xends(2)];
end

% Evaluate the function at the turning points:
y = feval(f, x);

end
