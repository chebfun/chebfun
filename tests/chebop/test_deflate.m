function pass = test_deflate(~)
% TEST_DEFLATE   Test the deflation capabilities of Chebfun

%% Bratu equation
N = chebop(@(x,u) diff(u,2) + 2*exp(u), [0 1], 0, 0);
r0 = N\0;
Ndef = deflate(N, r0, 1, 0);
r1 = Ndef\0;

% Ensure we got two solutions
pass(1, 1) = norm(N(r0)) < 1e-8 && norm(N(r1)) < 1e-8;
% Ensure we got different solutions
pass(2, 1) = norm(r0-r1) > 1;

%% Painleve
L = 10;
d = [0 L];
N = chebop(@(x,u) diff(u,2)-u.^2+x, d);

% Boundary conditions
N.lbc = 0;
N.rbc = sqrt(L);

r0 = N\0; % First solution
Ndef = deflate(N, r0, 3, .1); % Deflate for second solution;
r1 = Ndef\0;

% Ensure we got two solutions
pass(1, 2) = norm(N(r0)) < 1e-9 && norm(N(r1)) < 1e-9;
% Ensure we got different solutions
pass(2, 2) = norm(r0-r1) > 1;
%% Herceg
d = [0, 1];
ep = .25;
N = chebop(@(x,u) -ep^2*diff(u,2) + (u.^2 + u - 0.75).*(u.^2 + u - 3.75), d);
N.bc = 0;

x = chebfun('x', d);
u0 = 0*x;
N.init = u0;

r0 = N\0;
Ndef = deflate(N, r0, 1, 0);
r1 = Ndef\0;
% Deflate again
Ndef = deflate(N, [r0 r1], 1, 0);
r2 = Ndef\0;

% Ensure we got three solutions
pass(1, 3) = norm(N(r0)) < 1e-9 && norm(N(r1)) < 1e-9 && norm(N(r2)) < 1e-9;
% Ensure we got three different solutions
pass(2, 3) = norm(r0 - r1) > .1 && norm(r0 - r2) > .1 && norm(r1 - r2) > .1;
