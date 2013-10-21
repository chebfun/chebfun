function pass = test_definePoint(pref)

if ( nargin == 0 )
    pref = chebpref();
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

% Test an array-valued CHEBFUN object at a number of points:
f = chebfun(@(x) [x, x], [-1, 0, 1], pref);
% Use DEFINEPOINT:
f = definePoint(f, [-.25, .5], [1, 2 ; 3 4]);
pass(5) = all(size(f.impulses) == [5, 2]) && ...
    all(all(f.impulses == [-1 -1 ; 1 2 ; 0 0 ; 3 4 ; 1 1]));
% Use SUBSASGN:
f([-.25, .5]) = [1, 2 ; 3 4];
pass(6) = all(size(f.impulses) == [5, 2]) && ...
    all(all(f.impulses == [-1 -1 ; 1 2 ; 0 0 ; 3 4 ; 1 1]));

% Test an array-valued CHEBFUN object at a number of points (vector expansion):
f = chebfun(@(x) [x, x], [-1, 0, 1], pref);
% Use DEFINEPOINT:
f = definePoint(f, [-.25, .5], [1, 2]);
pass(7) = all(size(f.impulses) == [5, 2]) && ...
    all(all(f.impulses == [-1 -1 ; 1 2 ; 0 0 ; 1 2 ; 1 1]));
% Use SUBSASGN:
f([-.25, .5]) = [1, 2];
pass(8) = all(size(f.impulses) == [5, 2]) && ...
    all(all(f.impulses == [-1 -1 ; 1 2 ; 0 0 ; 1 2 ; 1 1]));

% Test an array-valued CHEBFUN object at a number of points (scalar expansion):
f = chebfun(@(x) [x, x], [-1, 0, 1], pref);
% Use DEFINEPOINT:
f = definePoint(f, [-.25, .5], 1);
pass(9) = all(size(f.impulses) == [5, 2]) && ...
    all(all(f.impulses == [-1 -1 ; 1 1 ; 0 0 ; 1 1 ; 1 1]));
% Use SUBSASGN:
f([-.25, .5]) = 1;
pass(10) = all(size(f.impulses) == [5, 2]) && ...
    all(all(f.impulses == [-1 -1 ; 1 1 ; 0 0 ; 1 1 ; 1 1]));

end
