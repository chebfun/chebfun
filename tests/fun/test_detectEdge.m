% Test file for @chebfun/detectEdge.m.

function pass = test_detectEdge(pref)

% We test detectEdge directly by giving it function handles with known
% discontinuities in the 0, ..., 4th derivatives.

if ( nargin == 0 )
    pref = chebfunpref;
end

% Initialise seed:
seedRNG(13);

% Number of points to test:
M = 10;
% Generate random points in [0, 1]:
x0 = rand(M,1);

% Initialise edge vector:
edge = zeros(M,1);

f = classicfun.constructor(0,struct('domain', [0 1]));

% Jump:
for j = 1:M
    edge(j) = fun.detectEdge(f, @(x) exp(x)+cos(7*x)+0.1*sign(x-x0(j)), 1, ...
        1, pref);
end
err = norm(edge - x0, inf);
pass(1) = err < 5e-14;

% C1:
for j = 1:M
    edge(j) = fun.detectEdge(f, @(x) exp(x)+cos(7*x)+0.1*abs(x-x0(j)), 1, ...
        1, pref);
end
err = norm(edge - x0, inf);
pass(2) = err < 5e-14;

% C2:
for j = 1:M
    edge(j) = fun.detectEdge(f, @(x) sign(x-x0(j)).*(x-x0(j)), 1, 1, pref);
end
err = norm(edge - x0, inf);
pass(3) = err < 5e-14;

% C3:
for j = 1:M
    edge(j) = fun.detectEdge(f, @(x) abs(x-x0(j)).^3, 1, 1, pref);
end
err = norm(edge - x0, inf);
pass(4) = err < 5e-14;

% C4:
for j = 1:M
    edge(j) = fun.detectEdge(f, @(x) sign(x-x0(j)).*(x-x0(j)).^3, 1, 1, pref);
end
err = norm(edge - x0, inf);
pass(5) = err < 5e-14;

%% Two examples from #809:
f = classicfun.constructor(0, struct('domain', [0 1]));
x0 = hex2num('3fe55ec001aabd80');  % 0.667816165214319
op = @(x) exp(x) + cos(7*x) + 0.1*sign(x - x0);
edge = fun.detectEdge(f, op, 1, 1, pref);
err = abs(edge - x0);
pass(6) = ( err < 5*eps );

% TODO:  #809 still causes this test to fail.  Comment it back in when that
% issue is resolved.
%f = classicfun.constructor(0, struct('domain', [0 1]));
%x0 = hex2num('3fce67c02dc03de0');
%op = @(x) exp(x)+cos(7*x)+abs(x-x0);
%edge = fun.detectEdge(f, op, 1, 1, pref);
%err = abs(edge - x0)
%pass(6) = ( err < 5*eps );

%% A test for a blow-up function with large domain. This is to exercise the
% relevant part in findBlowup@detectEdge.
f = classicfun.constructor(0, struct('domain', [1 2]));
op = @(x) tan(x);
pref.blowup = true;
edge = fun.detectEdge(f, op, 1, 1, pref);
err = abs(edge - pi/2);
pass(7) = ( err < 5*eps );

end
