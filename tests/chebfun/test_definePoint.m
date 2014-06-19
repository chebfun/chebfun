function pass = test_definePoint(pref)

if ( nargin == 0 )
    pref = chebfunpref();
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
pass(5) = all(size(f.pointValues) == [5, 2]) && ...
    all(all(f.pointValues == [-1 -1 ; 1 2 ; 0 0 ; 3 4 ; 1 1]));
% Use SUBSASGN:
f([-.25, .5]) = [1, 2 ; 3 4];
pass(6) = all(size(f.pointValues) == [5, 2]) && ...
    all(all(f.pointValues == [-1 -1 ; 1 2 ; 0 0 ; 3 4 ; 1 1]));

% Test an array-valued CHEBFUN object at a number of points (vector expansion):
f = chebfun(@(x) [x, x], [-1, 0, 1], pref);
% Use DEFINEPOINT:
f = definePoint(f, [-.25, .5], [1, 2]);
pass(7) = all(size(f.pointValues) == [5, 2]) && ...
    all(all(f.pointValues == [-1 -1 ; 1 2 ; 0 0 ; 1 2 ; 1 1]));
% Use SUBSASGN:
f([-.25, .5]) = [1, 2];
pass(8) = all(size(f.pointValues) == [5, 2]) && ...
    all(all(f.pointValues == [-1 -1 ; 1 2 ; 0 0 ; 1 2 ; 1 1]));

% Test an array-valued CHEBFUN object at a number of points (scalar expansion):
f = chebfun(@(x) [x, x], [-1, 0, 1], pref);
% Use DEFINEPOINT:
f = definePoint(f, [-.25, .5], 1);
pass(9) = all(size(f.pointValues) == [5, 2]) && ...
    all(all(f.pointValues == [-1 -1 ; 1 1 ; 0 0 ; 1 1 ; 1 1]));
% Use SUBSASGN:
f([-.25, .5]) = 1;
pass(10) = all(size(f.pointValues) == [5, 2]) && ...
    all(all(f.pointValues == [-1 -1 ; 1 1 ; 0 0 ; 1 1 ; 1 1]));

%% Test on singular function: piecewise smooth chebfun

% define the domain:
dom = [-2 -1 0 1];

op1 = @(x) sin(x);
op2 = @(x) 1./((1+x).^0.5);
op3 = @(x) x+1;
op = {op1, op2, op3};

f = chebfun(op, dom, 'exps', [0 0 -0.5 0 0 0]);
brkpts = zeros(1,2);
brkpts(1) = -0.7;
brkpts(2) = -0.5;
g = definePoint(f, brkpts(1), 1);
g(brkpts(2)) = 2;

% check values:
check = zeros(1,4);
check(1) = all(g.domain == unique([dom, brkpts]));
check(2) = feval(g, brkpts(1)) == 1;
check(3) = feval(g, brkpts(2)) == 2;
check(4) = all(g.pointValues == [f.pointValues(1:2); 1; 2; f.pointValues(3:4)]);

pass(11) = all( check );

%% Test for function defined on unbounded domain:

% Functions on [-inf inf]:

% Set the domain:
dom = [-Inf Inf];

% Blow-up function:
op = @(x) (1-exp(-x.^2));
f = chebfun(op, dom); 
brkpts = zeros(1,2);
brkpts(1) = -7;
brkpts(2) = 2;
g = definePoint(f, brkpts(1), 3);
g(brkpts(2)) = 4;

% check values:
check = zeros(1,4);
check(1) = all(g.domain == unique([dom, brkpts]));
check(2) = feval(g, brkpts(1)) == 3;
check(3) = feval(g, brkpts(2)) == 4;
check(4) = all(g.pointValues == [f.pointValues(1); 3; 4; f.pointValues(end)]);

pass(12) = all( check );

%% Test assigning to an endpoint (#896)
f = chebfun(0, [-1, 1]);
f(1) = 1;
pass(13) = abs(feval(f, 1) - 1) < epslevel(f);

end
