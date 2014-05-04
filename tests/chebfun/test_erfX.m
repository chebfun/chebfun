% Test file for all the functions related to the error function ERF().

function pass = test_erfX(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

F = {@erf, @erfinv, @erfc,  @erfcx, @(x) erfcinv(x+1)};

pass = false(5, 2);

for k = 1:numel(F)

    f = chebfun(@sin, pref);
    g = feval(F{k}, f);
    h = chebfun(@(x) feval(F{k}, sin(x)), pref);
    pass(k, 1) = normest(g - h) < 10*epslevel(h);
    
    f = chebfun(@(x) [sin(x), exp(x-1.1), x/2], pref);
    g = feval(F{k}, f);
    h = chebfun(@(x) [feval(F{k}, sin(x)), ...
                      feval(F{k}, exp(x-1.1)), ...
                      feval(F{k}, x/2)], pref);
    pass(k, 2) = normest(g - h) < 10*epslevel(h)*vscale(h);
    
end

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
