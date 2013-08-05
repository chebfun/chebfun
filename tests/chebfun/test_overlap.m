function pass = test_overlap(pref)

% Grab some preferences:
if ( nargin == 0 )
    pref = chebfun.pref();
end

f = chebfun(@sin, [-1 -.5 0 1], pref);
g = chebfun(@sin, [-1 0 0.5 1], pref);

[f2, g2] = overlap(f, g);
xx = linspace(-1, 1);
pass(1)  = norm(feval(f, xx) - feval(f2, xx), inf) < 10*f.epslevel.*f.vscale;
pass(2)  = norm(feval(g, xx) - feval(g2, xx), inf) < 10*g.epslevel.*g.vscale;

end
