% Test file for @chebfun/cummax.m.

function pass = test_cummax(pref)

% Check operation in the general case
f = chebfun('exp(x)*sin(10*x)');
g = cummax(f);
pass(1) = abs(min(g-f)) < 10*vscale(f)*eps;
pass(2) = abs(max(g)-max(f)) < 10*vscale(f)*eps;

% Check the output on an increasing function
f = chebfun('exp(x)');
g = cummax(f);
pass(3) = abs(max(f-g)) < 10*vscale(f)*eps;

% Check the behavior on a row chebfun
f = f';
g = cummax(f);
pass(4) = abs(max(f-g)) < 10*vscale(f)*eps;

% Check error condition
x = chebfun('x');
f = [x x.^2];
try
    g = cummax(f);
    pass(5) = false;
catch ME
    pass(5) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:cummax:quasi');

end