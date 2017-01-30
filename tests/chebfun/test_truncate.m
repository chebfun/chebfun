function pass = test_truncate(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

%% Test truncate():
x = chebfun('x', pref);
f = sign(x);
g = truncate(f, 10);
pass(1) = length(g) == 10;
c = chebcoeffs(g);
pass(2) = abs(c(end) - 4/(9*pi)) < 10*eps;

%% Test 'trunc' flag:
G = chebfun(@(x) sign(x), 'trunc', 10, pref);
pass(3) = length(g) == 10;
pass(4) = norm(g - G, inf) < 10*eps;
c = chebcoeffs(g);
pass(5) = abs(c(end) - 4/(9*pi)) < 10*eps;

%% Test trig-based:
f = chebfun(@(x) exp(sin(pi*x)), 'trig', pref);
g = truncate(f, 10);
pass(6) = length(g) == 10;
tech = get(g.funs{1}, 'tech');
pass(7) = isa(tech(), 'trigtech');

%% Test trig-based on a non-standard interval:
f = chebfun('exp(sin(t))', [0, 2*pi]);
p = chebfun(f,'trunc', 31,'trig');
pass(8) = all(domain(p) == [0, 2*pi]);
pass(9) = norm(f-p) < 100*eps;

end
