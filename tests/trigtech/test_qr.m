% Test file for trigtech/qr.m

function pass = test_qr(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

testclass = trigtech();

%%
% Do a few spot-checks.

f = testclass.make(@(x) exp(sin(pi*x)), [], pref);
pass(1:2) = test_one_qr(f, x);
pass(3:4) = test_one_qr_with_perm(f, x);

f = testclass.make(@(x) [exp(sin(pi*x)) 3./(4-cos(pi*x))], [], pref);
pass(5:6) = test_one_qr(f, x);
pass(7:8) = test_one_qr_with_perm(f, x);

f = testclass.make(@(x) [ones(size(x)) sin(pi*x) cos(pi*x) sin(2*pi*x) ...
                         cos(2*pi*x) sin(3*pi*x)], [], pref);
pass(9:10) = test_one_qr(f, x);
pass(11:12) = test_one_qr_with_perm(f, x);

f = testclass.make(@(x) [3./(4-exp(1i*pi*x)) exp(sin(pi*x)) cos(3*pi*x)], ...
[], pref);
pass(13:14) = test_one_qr(f, x);
pass(15:16) = test_one_qr_with_perm(f, x);

%%
% Check that the 'vector' flag works properly.
N = size(f, 2);
[Q1, R1, E1] = qr(f, []);
[Q2, R2, E2] = qr(f, 'vector');
err = E1(:, E2) - eye(N);
pass(17) = all(err(:) == 0);



%%
% Check a rank-deficient problem:
% [TODO]: Is this correct?
%
% Rank deficient QR factorizations fail due to the bug 
% in issue #1441. These tests are disabled until the bug
% is addressed.
f = testclass.make(@(x) [cos(pi*x) cos(pi*x) cos(pi*x)], [], pref);
[Q, R] = qr(f, []);
Q = simplify(Q,100*eps);
%pass(18) = all(size(Q) == 3) && all(size(R) == 3);
pass(18) = 1;
I = eye(3);
%pass(19) = norm(innerProduct(Q, Q) - I, inf) < ...
%10*max(vscale(f)*eps);
pass(19) = 1;
% These tests should be reverted once issue #1441 is
% fixed.



%%
% Check that the vscale comes out with the correct size for QR of an
% array-valued chebtech.
f = testclass.make(@(x) [sin(pi*x) cos(pi*x) cos(2*pi*x)], [], pref);
[Q, R] = qr(f, []);
pass(20) = isequal(size(vscale(Q)), [1 3]);

end

% Tests the QR decomposition for a CHEBTECH object F using a grid of points X
% in [-1  1] for testing samples.
function result = test_one_qr(f, x)
    N = size(f, 2);
    [Q, R] = qr(f);

    % Check orthogonality.
    ip = innerProduct(Q, Q);
    result(1) = max(max(abs(ip - eye(N)))) < 10*max(vscale(f)*eps);

    % Check that the factorization is accurate.
    err = Q*R - f;
    result(2) = norm(feval(err, x), inf) < 100*max(vscale(f)*eps);
end

% Same as the previous function but this time uses the QR factorization with
% permutations.
function result = test_one_qr_with_perm(f, x)
    N = size(f, 2);
    [Q, R, E] = qr(f);

    % Check orthogonality.
    ip = innerProduct(Q, Q);
    result(1) = max(max(abs(ip - eye(N)))) < 10*max(vscale(f)*eps);

    % Check that the factorization is accurate.
    err = Q*R - f*E;
    result(2) = norm(feval(err, x), inf) < 100*max(vscale(f)*eps);
end
