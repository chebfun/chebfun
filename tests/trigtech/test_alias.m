% Test file for trigtech/alias.m

function pass = test_alias(varargin)

% Set a tolerance (pref.eps doesn't matter)
tol = 10*eps;

%% 
% Test padding of a vector of coefficients.

c0 = (1:9)';
c1 = trigtech.alias(c0, 13);
pass(1) = norm([0; 0; c0; 0; 0] - c1, inf) == 0;

%% 
% Compare against result of evaluating on a smaller grid:
n = [1 2 3 7 10 12 15];
start = 1;
for j=1:numel(n)
    pass(start+j) = test_alias_by_interpolating(@(x) cos(1+sin(2*pi*x)), n(j), tol);
end

% Complex function:
n = [1 2 3 7 10 12 15];
start = length(pass);
for j=1:numel(n)
    pass(start+j) = test_alias_by_interpolating(@(x) cos(1+sin(pi*x)) + 1i*exp(cos(pi*x)), n(j), tol);
end

% Larger grid:
n = [8 15];
start = length(pass);
for j=1:numel(n)
    pass(start+j) = test_alias_by_interpolating(@(x) 1 + cos(3*pi*x), n(j), tol);
end

% Complex function:
n = [8 15];
start = length(pass);
for j=1:numel(n)
    pass(start+j) = test_alias_by_interpolating(@(x) 1 + cos(3*pi*x) + 1i*sin(2*pi*x), n(j), tol);
end

% Array-valued function:
start = length(pass);
n = [1 2 3 7 10 12 15];
for j=1:numel(n)
    pass(start+j) = test_alias_by_interpolating(@(x) [cos(pi*x) sin(1+cos(pi*x)) exp(cos(pi*x))], n(j), tol);    
end

end

function pass = test_alias_by_interpolating(op, n, tol)

f = trigtech.make(op);
c = f.alias(f.coeffs, n);
x = trigtech.trigpts(size(c, 1));
cexact = f.vals2coeffs(op(x));

pass = norm(c - cexact, inf) < tol;

end
