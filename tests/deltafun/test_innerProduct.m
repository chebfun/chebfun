% Test file for @deltafun/sum.m

function pass = test_innerProduct(pref)

if (nargin < 1)
    pref = chebfunpref();
end

% Get the tolerance:
tol = pref.deltaPrefs.deltaTol;
d = deltafun();

f = deltafun(fun.constructor(@(x) exp(x)), []);
zf = fun.constructor( @(x) 0.*x );
mag = rand(1,5);
loc = rand(1,5);

g = deltafun(f.funPart, struct('deltaMag', mag, 'deltaLoc', loc));

pass(1) = isempty(innerProduct(d, d)) && isempty(innerProduct(d, g)) && ...
          isempty(innerProduct(g, d));

pass(2) = abs(innerProduct(g, f) - (sum(mag.*feval(f, loc)) + sum(f.*f))) < tol;

d = deltafun(zf, struct('deltaMag', [0; 1], 'deltaLoc', 0) );
g = deltafun(f);

pass(3) = abs(innerProduct(d, g) + feval(diff(g), 0)) < tol;

pass(4) = abs(innerProduct(diff(d), g) - feval(diff(diff(g)), 0)) < tol;

loc = [-.5, 0, .5];
mag = [1 2 3; -1 0 1; 0 0 1];
d = deltafun(zf, struct('deltaMag', mag, 'deltaLoc', loc));
ip = sum(mag(1,:).*feval(f, loc)) - sum(mag(2,:).*feval(diff(f), loc)) + ...
    sum(mag(3,:).*feval(diff(f,2), loc));
pass(5) = abs(innerProduct(g, d) - ip) < tol;

d = deltafun(zf, struct('deltaMag', 1, 'deltaLoc', 0));
ip = innerProduct(d, d);
pass(6) = isinf(ip) && ip > 0;

ip = innerProduct(d, -d);
pass(7) = isinf(ip) && ip < 0;

end
