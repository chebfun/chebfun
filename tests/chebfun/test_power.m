function pass = test_power(pref)

if ( nargin == 0 ) 
    pref = chebpref();
end

% [TODO]: So far, POWER only supports easy cases (i.e., positive integers).
% Tests will have to expanded once this functionality is available.


%% Scalar-valued
f = chebfun(@(x) sin(x));
g = f.^0;
pass(1) = normest(g - 1) < 10*epslevel(f);

g = f.^1;
pass(2) = normest(g - f) < 10*epslevel(f);

g = f.^2;
h = chebfun(@(x) sin(x).^2);
pass(3) = normest(g - h) < 10*epslevel(h);

g = f.^3;
h = chebfun(@(x) sin(x).^3);
pass(4) = normest(g - h) < 10*epslevel(h);

%% Array-valued
f = chebfun(@(x) [sin(x), cos(x), 1i*exp(x)]);
g = f.^0;
pass(5) = min(size(g)) == 3 && normest(g - 1) < epslevel(f);

g = f.^1;
pass(6) = min(size(g)) == 3 && normest(g - f) < epslevel(f);

g = f.^2;
h = chebfun(@(x) [sin(x).^2, cos(x).^2, -exp(2*x)]);
pass(7) = min(size(g)) == 3 && normest(g - h) < 10*vscale(h)*epslevel(h);

g = f.^3;
h = chebfun(@(x) [sin(x).^3, cos(x).^3, -1i*exp(3*x)]);
pass(8) = min(size(g)) == 3 && normest(g - h) < 10*vscale(h)*epslevel(h);

%% constant .^ CHEBFUN
f = chebfun(@(x) sin(x));
g = 1.^f;
pass(9) = normest(g - 1) < 10*epslevel(f);

f = chebfun(@(x) sin(x));
g = (2i).^f;
h = chebfun(@(x) (2i).^sin(x));
pass(10) = normest(g - h) < 10*epslevel(f);

f = chebfun(@(x) [sin(x), cos(x), 1i*exp(x)]);
g = 1.^f;
pass(11) = min(size(g)) == 3 && normest(g - 1) < 10*epslevel(h);

g = (2i).^f;
h = chebfun(@(x) [2i.^sin(x), 2i.^cos(x), 2i.^(1i*exp(x))]);
pass(12) = min(size(g)) == 3 && normest(g - h) < 100*epslevel(h);

%% CHEBFUN .^ CHEBFUN

x = chebfun(@(x) x, [.1, 2]);
g = x.^x;
h = chebfun(@(x) x.^x, [.1, 2]);
pass(13) = normest(g - h, [.1, 2]) < 10*epslevel(h)*vscale(h);

x = chebfun(@(x) [x, exp(1i*x)], [.1, 2]);
g = x.^x;
h = chebfun(@(x) [x.^x, exp(1i*x).^exp(1i*x)], [.1, 2]);
pass(14) = normest(g - h, [.1, 2]) < 10*epslevel(h)*vscale(h);

end

function out = normest(f, dom)

% Generate a few random points to use as test values.
seedRNG(6178);
if ( nargin == 1 )
    x = 2 * rand(100, 1) - 1;
else
    x = sum(dom) * rand(10, 1) - dom(1);
end

out = norm(feval(f, x), inf);

end

