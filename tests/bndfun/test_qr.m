% Test file for bndfun/qr.m

function pass = test_qr(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% Set the domain
dom = [-2 7];

% Create seed for random number generator
seedRNG(6178);

% Generate a few random points to use as test values.
x = diff(dom) * rand(100, 1) + dom(1);

%%
% Do a few spot-checks.

f = bndfun(@(x) sin(x), struct('domain', dom), pref);
pass(1:2) = test_one_qr(f, x);
pass(3:4) = test_one_qr_with_perm(f, x);

f = bndfun(@(x) [cos(x) exp(x)], struct('domain', dom), pref);
pass(5:6) = test_one_qr(f, x);
pass(7:8) = test_one_qr_with_perm(f, x);

f = bndfun(@(x) [ones(size(x)) x x.^2 x.^3 x.^4 x.^5 x.^6 x.^7], ...
    struct('domain', dom), pref);
pass(9:10) = test_one_qr(f, x);
pass(11:12) = test_one_qr_with_perm(f, x);

f = bndfun(@(x) [1./(1+1i*x.^2) sinh((1-1i)*x) (exp(x) - x.^3)], ...
    struct('domain', dom), pref);
pass(13:14) = test_one_qr(f, x);
pass(15:16) = test_one_qr_with_perm(f, x);

%%
% Check that the 'vector' flag works properly.
N = size(f, 2);
[Q1, R1, E1] = qr(f);
[Q2, R2, E2] = qr(f, 'vector');
err = E1(:, E2) - eye(N);
pass(17) = all(err(:) == 0);

%%
% Check a rank-deficient problem:
f = bndfun(@(x) [x x x], struct('domain', dom), pref);
[Q, R] = qr(f);
pass(18) = all(size(Q) == 3) && all(size(R) == 3);
I = eye(3);
pass(19) = norm(innerProduct(Q, Q) - I, inf) < ...
    max(get(f, 'vscale').*get(f, 'epslevel'));

end

% Tests the QR decomposition for a BNDFUN object F using a grid of points X
% in [-1  1] for testing samples.
function result = test_one_qr(f, x)
    N = size(f, 2);
    [Q, R] = qr(f);

    % Check orthogonality.
    ip = innerProduct(Q, Q);
    result(1) = max(max(abs(ip - eye(N)))) < ...
        10*max(get(f, 'vscale').*get(f, 'epslevel'));

    % Check that the factorization is accurate.
    err = Q*R - f;
    result(2) = norm(feval(err, x), inf) < ...
        10*max(get(f, 'vscale').*get(f, 'epslevel'));
end

% Same as the previous function but this time uses the QR factorization with
% permutations.
function result = test_one_qr_with_perm(f, x)
    N = size(f, 2);
    [Q, R, E] = qr(f);

    % Check orthogonality.
    ip = innerProduct(Q, Q);
    result(1) = max(max(abs(ip - eye(N)))) < ...
        10*max(get(f, 'vscale').*get(f, 'epslevel'));

    % Check that the factorization is accurate.
    err = Q*R - f*E;
    result(2) = norm(feval(err, x), inf) < ...
        10*max(get(f, 'vscale').*get(f, 'epslevel'));
end
