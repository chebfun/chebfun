% Test file for @deltafun/sum.m

function pass = test_innerProduct(pref)

if (nargin < 1)
    pref = chebfunpref();
end
%%
% Get the tolerance:
tol = 1e-12;

d = deltafun();

f = fun.constructor(@(x) exp(x));
zf = fun.constructor( @(x) 0.*x );
mag = rand(1,5);
loc = rand(1,5);

g = deltafun(f, mag, loc);

pass(1) = isempty(innerProduct(d, d));

pass(2) = isempty(innerProduct(d, g));

pass(3) = abs(innerProduct(g, f) - (sum(mag.*feval(f, loc)) + sum(f.*f))) < tol;


d = deltafun(zf, [0; 1], 0 );
g = deltafun(f, [], [] );

pass(4) = abs(innerProduct(d, g) + feval(diff(g), 0)) < tol;

pass(5) = abs(innerProduct(diff(d), g) - feval(diff(diff(g)), 0)) < tol;

loc = [-.5, 0, .5];
mag = [1 2 3; -1 0 1; 0 0 1];
d = deltafun(zf, mag, loc);
ip = sum(mag(1,:).*feval(f, loc)) - sum(mag(2,:).*feval(diff(f), loc)) + ...
    sum(mag(3,:).*feval(diff(f,2), loc));
pass(6) = abs(innerProduct(g, d) - ip) < tol;

end