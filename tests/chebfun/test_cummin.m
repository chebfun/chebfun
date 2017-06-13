% Test file for @chebfun/cummin.m.

function pass = test_cummin(pref)

% Check operation in the general case
f = chebfun('exp(x)*sin(10*x)');
g = cummin(f);
pass(1) = abs(min(f-g)) < 10*vscale(f)*eps;
pass(2) = abs(min(f)-min(g)) < 10*vscale(f)*eps;

% Check the output on a decreasing function
f = chebfun('exp(1./(x+2))');
g = cummin(f);
pass(3) = abs(max(f-g)) < 10*vscale(f)*eps;

% Check the behavior on a row chebfun
f = f';
g = cummin(f);
pass(4) = abs(max(f-g)) < 10*vscale(f)*eps;

% Check error condition
x = chebfun('x');
f = [x x.^2];
try
    g = cummin(f);
    pass(5) = false;
catch ME
    pass(5) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:cummin:quasi');

end

end