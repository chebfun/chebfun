function pass = test_deflate_painleve(~)
% TEST_DEFLATE   Test the deflation capabilities of Chebfun

% Painleve

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
pass(1,1) = norm(N(r0)) < 1e-9 && norm(N(r1)) < 1e-9;
% Ensure we got different solutions
pass(1,2) = norm(r0-r1) > 1;

end
