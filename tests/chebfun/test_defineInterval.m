function pass = test_defineInterval(pref)

if ( nargin == 0 )
    pref = chebpref();
end

% Test a scalar-valued CHEBFUN object:
f = chebfun(@(x) x, [-1, 0, 1], pref);
% Use DEFINEINTERVAL:
f = defineInterval(f, [0 .5], 1);
pass(1) = all(f.domain == [-1, 0, .5, 1]) && feval(f, .25) == 1;
% Use SUBSASGN:
f{-.5, 0} = 2;
pass(2) = all(f.domain == [-1, -.5, 0, .5, 1]) && feval(f, -.25) == 2;

% Test an array-valued CHEBFUN object:
f = chebfun(@(x) [x, x], [-1, 0, 1], pref);
% Use DEFINEINTERVAL:
f = defineInterval(f, [0, .5], [1, 2]);
pass(3) = all(f.domain == [-1, 0, .5, 1]) && all(feval(f, .25) == [1, 2]);
% Use SUBSASGN:
f{-.5,0} = [2, 3];
pass(4) = all(f.domain == [-1, -.5, 0, .5, 1]) && all(feval(f, -.25) == [2, 3]);

% Test an array-valued CHEBFUN object (scalar expansion):
f = chebfun(@(x) [x, x], [-1, 1], pref);
% Use DEFINEINTERVAL:
f = defineInterval(f, [0, .5], 1);
pass(5) = all(f.domain == [-1, 0, .5, 1]) && all(feval(f, .25) == 1);
% Use SUBSASGN:
f{-.5,0} = -1;
pass(6) = all(f.domain == [-1, -.5, 0, .5, 1]) && all(feval(f, -.25) == -1);

% Test an array-valued CHEBFUN object (additional breaks):
f = chebfun(@(x) [x, x], [-1, 1], pref);
% Use DEFINEINTERVAL:
f = defineInterval(f, [0, .25, .5], 1);
pass(7) = all(f.domain == [-1, 0, .25, .5, 1]) && all(feval(f, .125) == 1);
% Use SUBSASGN:
f{-.5,-.1,0} = -1;
pass(8) = all(f.domain == [-1, -.5, -.1, 0, .25, .5, 1]) && ...
    all(feval(f, -.125) == -1);

% Test an scalar-valued CHEBFUN object redefined by a CHEBFUN.
f = chebfun(@(x) x, [-1, 1], pref);
g = chebfun(@(x) sin(x), [-.5, 0, .75], pref);
% Use DEFINEINTERVAL:
f = defineInterval(f, [-.5, .25], g);
pass(9) = all(f.domain == [-1, -.5, 0, .25, 1]) && ...
    norm(feval(f, .125) - sin(.125)) < 10*epslevel(f) && ...
    norm(feval(f, .625) - .625) < 10*epslevel(f);
% Use SUBSASGN:
f{-.5, .25} = g;
pass(10) = all(f.domain == [-1, -.5, 0, .25, 1]) && ...
    norm(feval(f, .125) - sin(.125)) < 10*epslevel(f) && ...
    norm(feval(f, .625) - .625) < 10*epslevel(f);

% Test an array-valued CHEBFUN object redefined by a CHEBFUN.
f = chebfun(@(x) [x, -x], [-1, 0, 1], pref);
g = chebfun(@(x) [-2*x, 2*x], [-.5, .5], pref);
% Use DEFINEINTERVAL:
f = defineInterval(f, [-.5, .5], g);
pass(11) = norm(feval(f, .25) - [-.5, .5]) < 10*epslevel(f);
% Use SUBSASGN:
f{-.5, .5} = g;
pass(12) = norm(feval(f, .25) - [-.5, .5]) < 10*epslevel(f);

% Test removal:
f = chebfun(@(x) x, pref);
f = defineInterval(f, [-.5, .5], []);
pass(13) = all(f.domain == [-1, -.5, 0]);
f = chebfun(@(x) x, pref);
f = defineInterval(f, [.5, 1], []);
pass(14) = all(f.domain == [-1, .5]);
f = chebfun(@(x) x, pref);
f = defineInterval(f, [-1, .25], []);
pass(15) = all(f.domain == [.25, 1]);

end
