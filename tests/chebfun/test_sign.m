function pass = test_sign(pref)

if ( nargin == 0 )
    pref = chebfun.pref();
end

% Pre-allocate pass matrix:
pass = zeros(3, 4); 

% Initialise random vector:
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

%% Simple tests
pref.chebfun.splitting = 0;
f = chebfun('x.^2', pref);   
tol = get(f, 'epslevel')*get(f, 'hscale');
f1 = sign(f);
pass(1,1) = all(feval(f1, x) == 1);

% Test if impulses are dealt with correctly: 
f = chebfun(@(x) cos(pi*x) + 2, sort(x)', pref);
f.impulses(:,:,1) = -pi;
f1 = sign(f);
pass(1,2) = all(f1.impulses(:,:,1) == -1);

% Test also on longer intervals:
f = chebfun(@(x) x, [-1, 100], pref);
f1 = sign(f);
h = chebfun({-1, 1}, [-1, 0, 100]);
pass(1,3) = normest(h - f1) < tol;

% Test functions with breakpoints:
f = chebfun(@(x) x, [-1, 1e-10, 1], pref);
f1 = sign(f);
h = restrict(h,[-1, 1]);
pass(1,4) = normest(f1 - h) < tol;

% Real, imaginary and complex CHEBFUN objects:
f = chebfun({-1, 1, -1, 1, -1, 1}, -3:3);

%% Real CHEBFUN
gHandle1 = @(x) sin(pi*x);
g = chebfun(@(x) gHandle1(x), -3:3, pref);
g1 = sign(g);
pass(2,1) = length(g1.funs) == 6;
pass(2,2) = normest(f - g1) < tol;
pref.chebfun.splitting = 1;
h1 = chebfun(@(x) sign( gHandle1(x) ), -3:3, pref);
pass(2,3) = length(h1.funs) == 6;
pass(2,4) = normest(f - h1) < tol;

%% Test complex CHEBFUN once that is supported?

%% Array-valued CHEBFUN
f1 = chebfun(@(x) feval(f, [x, x]) , -3:3, pref);
gHandle2 = @(x) cos(pi*(x-.5));
g = chebfun(@(x) [gHandle1(x), gHandle2(x)] , -3:3, pref);
g4 = sign(g);
pass(3,1) = length(g4.funs) == 6;
pass(3,2) = normest(f1 - g4) < tol;
h4 = chebfun(@(x) sign( [gHandle1(x), gHandle2(x)] ), -3:3, pref);
pass(3,3) = length(h4.funs) == 6;
pass(3,4) = normest(f1 - h4) < tol;

end