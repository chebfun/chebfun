% Test file for @chebfun/overlap.m.

function pass = test_overlap(pref)

% Grab some preferences:
if ( nargin == 0 )
    pref = chebpref();
end

% Test empty input.
f = chebfun();
g = chebfun();
[f2, g2] = overlap(f, g);
pass(1) = isempty(f2) && isempty(g2);

% Check behavior for input chebfuns with differing domains.
f = chebfun(@sin, [-1 -.5 0 1], pref);
g = chebfun(@sin, [-2 0 0.5 2], pref);

try
    [f2, g2] = overlap(f, g);
    pass(2) = false;
catch ME
    pass(2) = strcmp(ME.identifier, 'CHEBFUN:overlap:domains');
end

% Check behavior in the basic case.
f = chebfun(@sin, [-1 -.5 0 0.5 1], pref);
g = chebfun(@sin, [-1 0 0.5 1], pref);

[f2, g2] = overlap(f, g);
xx = linspace(-1, 1);
pass(3) = isequal(f2.domain, g2.domain) && ...
    norm(feval(f, xx) - feval(f2, xx), inf) < 10*epslevel(f)*vscale(f) && ...
    norm(feval(g, xx) - feval(g2, xx), inf) < 10*epslevel(g)*vscale(g);

% [TODO]: Are these the same tests as above?
% Check correct behavior for higher-order impulses.  [TODO]:  Use a function
% with real higher-order impulses instead of just faking them like this.
[f2, g2] = overlap(f, g);
xx = linspace(-1, 1);
pass(4) = isequal(f2.domain, g2.domain) && ...
    norm(feval(f, xx) - feval(f2, xx), inf) < 10*epslevel(f)*vscale(f) && ...
    norm(feval(g, xx) - feval(g2, xx), inf) < 10*epslevel(g)*vscale(g);

[g2, f2] = overlap(g, f);
xx = linspace(-1, 1);
pass(5) = isequal(f2.domain, g2.domain) && ...
    norm(feval(f, xx) - feval(f2, xx), inf) < 10*epslevel(f)*vscale(f) && ...
    norm(feval(g, xx) - feval(g2, xx), inf) < 10*epslevel(g)*vscale(g);
end
