function pass = test_mean(pref)

if ( nargin == 0 )
    pref = chebfun.pref();
end

f = chebfun(@sin, pref);
g = chebfun(@cos, pref);
h = chebfun(@(x) .5*(sin(x) + cos(x)), pref);
pass(1) = normest(mean(f, g) - h) < 10*epslevel(h);

f = chebfun(@(x) [sin(x), cos(x)], pref);
g = chebfun(@(x) [cos(x), 1i*exp(x)], pref);
h = .5*(f+g);
pass(2) = normest(mean(f, g) - h) < 10*epslevel(h);

% [TODO]: This requires SUM() before tests can be written
% f = chebfun(@sin, pref);
% pass(3) = abs(mean(f)) < epslevel(f);
pass(3) = false;

% [TODO]: More tests of this type.

% [TODO]: Test unbounded domains.

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