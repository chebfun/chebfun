% Test file for @chebfun/find.m.

function pass = test_find(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Check the empty case:
pass(1) = isempty(find(chebfun()));

% Check a few simple examples:
f = chebfun(@sin, [0 2*pi], pref);
x = find(f == 1/2) / pi;
pass(2) = norm(x - [1/6 ; 5/6], inf) < 10*vscale(f)*epslevel(f);

f = chebfun(@exp, [-1 -0.5 0 0.5 1], pref);
x = find(f == 1);
pass(3) = (length(x) == 1) && (abs(x) < 10*vscale(f)*epslevel(f));
x = find(f <= 0);
pass(4) = isempty(x);

f = chebfun(@(x) [sin(x) exp(x)], [-pi pi], pref);
[x, col] = find(~logical(f - 0.5)); % (f == 0.5 won't work for array-valued f.)
pass(5) = (length(x) == 3) && (length(col) == 3);
fx = feval(f, x);
err = [fx(1, col(1)) ; fx(2, col(2)) ; fx(3, col(3))];
pass(6) = norm(err - 0.5, inf) < 10*vscale(f)*epslevel(f);

% Check behavior for row chebfuns.
[row, y] = find(~logical(f.' - 0.5));
pass(7) = isequal(row, col) && isequal(x, y);

% Check error conditions.
try
    x = find(f);
    pass(8) = false;
catch ME
    pass(8) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:find:arrout');
end

try
    f = chebfun(@sin, pref);
    g = find(f);
    pass(9) = false;
catch ME
    pass(9) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:find:infset');
end

end
