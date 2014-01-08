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
fq = quasimatrix(f);

g = f.^0;
pass(5) = min(size(g)) == 3 && normest(g - 1) < epslevel(f);
gq = fq.^0;
pass(6) = normest(gq - g) < epslevel(g)*vscale(g);

g = f.^1;
pass(7) = min(size(g)) == 3 && normest(g - f) < epslevel(f);
gq = fq.^1;
pass(8) = normest(gq - g) < epslevel(g)*vscale(g);

g = f.^2;
h = chebfun(@(x) [sin(x).^2, cos(x).^2, -exp(2*x)]);
pass(9) = min(size(g)) == 3 && normest(g - h) < 10*vscale(h)*epslevel(h);
gq = fq.^2;
pass(10) = normest(gq - g) < epslevel(g);

g = f.^3;
h = chebfun(@(x) [sin(x).^3, cos(x).^3, -1i*exp(3*x)]);
pass(11) = min(size(g)) == 3 && normest(g - h) < 10*vscale(h)*epslevel(h);
gq = fq.^3;
pass(12) = normest(gq - g) < epslevel(g)*vscale(g);

%% constant .^ CHEBFUN
f = chebfun(@(x) sin(x));
g = 1.^f;
pass(13) = normest(g - 1) < 10*epslevel(f);

f = chebfun(@(x) sin(x));
g = (2i).^f;
h = chebfun(@(x) (2i).^sin(x));
pass(14) = normest(g - h) < 10*epslevel(f);

f = chebfun(@(x) [sin(x), cos(x), 1i*exp(x)]);
fq = quasimatrix(f);
g = 1.^f;
pass(15) = min(size(g)) == 3 && normest(g - 1) < 10*epslevel(h);
gq = 1.^fq;
pass(16) = normest(gq - g) < 10*epslevel(h)*vscale(h);;

g = (2i).^f;
h = chebfun(@(x) [2i.^sin(x), 2i.^cos(x), 2i.^(1i*exp(x))]);
pass(17) = min(size(g)) == 3 && normest(g - h) < 100*epslevel(h);
gq = (2i).^fq;
pass(18) = normest(gq - g) < 10*epslevel(h)*vscale(h);;

%% CHEBFUN .^ CHEBFUN

x = chebfun(@(x) x, [.1, 2]);
g = x.^x;
h = chebfun(@(x) x.^x, [.1, 2]);
pass(19) = normest(g - h, [.1, 2]) < 10*epslevel(h)*vscale(h);

x = chebfun(@(x) [x, exp(1i*x)], [.1, 2]);
g = x.^x;
h = chebfun(@(x) [x.^x, exp(1i*x).^exp(1i*x)], [.1, 2]);
pass(20) = normest(g - h, [.1, 2]) < 10*epslevel(h)*vscale(h);

xq = quasimatrix(x);
gq = xq.^xq;
pass(21) = normest(g - gq, [.1, 2]) < 10*epslevel(h)*vscale(h);

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

