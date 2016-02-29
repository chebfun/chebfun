function pass = test_deflate_bratu(~)
% TEST_DEFLATE   Test the deflation capabilities of Chebfun

% Bratu equation
N = chebop(@(x,u) diff(u,2) + 2*exp(u), [0 1], 0, 0);
r0 = N\0;
Ndef = deflate(N, r0, 1, 0);
r1 = Ndef\0;

% Ensure we got two solutions
pass(1,1) = norm(N(r0)) < 1e-8 && norm(N(r1)) < 1e-8;
% Ensure we got different solutions
pass(1,2) = norm(r0-r1) > 1;

end