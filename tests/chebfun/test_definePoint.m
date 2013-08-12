function pass = test_definePoint(pref)

if ( nargin == 0 )
    pref = chebfun.pref();
end

% Test a scalar-valued CHEBFUN object:
f = chebfun(@(x) x, [-1, 0, 1], pref);
% Use DEFINEPOINT:
f = definePoint(f, .5, 1);
pass(1) = all(f.domain == [-1, 0, .5, 1]) && feval(f, .5) == 1;
% Use SUBSASGN:
f(-.5) = 2;
pass(2) = all(f.domain == [-1, -.5, 0, .5, 1]) && feval(f, -.5) == 2;

% Test an array-valued CHEBFUN object:
f = chebfun(@(x) [x, x], [-1, 0, 1], pref);
% Use DEFINEPOINT:
f = definePoint(f, .5, [1, 2]);
pass(3) = all(f.domain == [-1, 0, .5, 1]) && all(feval(f, .5) == [1, 2]);
% Use SUBSASGN:
f(-.5) = [2, 3];
pass(4) = all(f.domain == [-1, -.5, 0, .5, 1]) && all(feval(f, -.5) == [2, 3]);

end
