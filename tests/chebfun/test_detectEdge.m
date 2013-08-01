function pass = test_detectEdge(pref)

% We test detectEdge directly by giving it function handles with known
% discontinuities in the 0, ..., 4th derivatives.

% Initialise seed:
rng(13);

% Number of points to test:
M = 10;
% Generate random points in [0, 1]:
x0 = rand(M,1);

% Initialise edge vector:
edge = zeros(M,1);

% Jump:
for j = 1:M
    edge(j) = chebfun.detectEdge(@(x) exp(x)+cos(7*x)+0.1*sign(x-x0(j)), [0, 1]);
end
pass(1) = norm(edge - x0, inf) < 5e14;

% C1:
for j = 1:M
    edge(j) = chebfun.detectEdge(@(x) exp(x)+cos(7*x)+0.1*abs(x-x0(j)), [0, 1]);
end
pass(2) = norm(edge - x0, inf) < 5e14;

% C2:
for j = 1:M
    edge(j) = chebfun.detectEdge(@(x) sign(x-x0(j)).*x.^2, [0, 1]);
end
pass(3) = norm(edge - x0, inf) < 5e14;

% C3:
for j = 1:M
    edge(j) = chebfun.detectEdge(@(x) abs(x-x0(j)).^3, [0, 1]);
end
pass(4) = norm(edge - x0, inf) < 5e14;

% C4:
for j = 1:M
    edge(j) = chebfun.detectEdge(@(x) sign(x-x0(j)).*x.^4, [0, 1]);
end
pass(5) = norm(edge - x0, inf) < 5e14;

end

