function pass = test_deflate_herceg(~)
% TEST_DEFLATE   Test the deflation capabilities of Chebfun

%% herceg equation
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
pass(1,1) = norm(N(r0)) < 1e-9 && norm(N(r1)) < 1e-9 && norm(N(r2)) < 1e-9;
% Ensure we got three different solutions
pass(1,2) = norm(r0 - r1) > .1 && norm(r0 - r2) > .1 && norm(r1 - r2) > .1;

end
