% Test file for chebtech/qr.m

function pass = test_qr(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebtech.pref;
end

% Generate a few random points to use as test values.
rngstate = rng();
rng(6178);
x = 2 * rand(100, 1) - 1;

for ( n = 1:2 )
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    %%
    % Do a few spot-checks.
    
    f = testclass.make(@(x) sin(x), [], [], pref);
    pass(n, 1:2) = test_one_qr(f, x);
    pass(n, 3:4) = test_one_qr_with_perm(f, x);
    
    f = testclass.make(@(x) [cos(x) exp(x)], [], [], pref);
    pass(n, 5:6) = test_one_qr(f, x);
    pass(n, 7:8) = test_one_qr_with_perm(f, x);
    
    f = testclass.make(@(x) [ones(size(x)) x x.^2 x.^3 x.^4 x.^5 x.^6 x.^7], ...
        [], [], pref);
    pass(n, 9:10) = test_one_qr(f, x);
    pass(n, 11:12) = test_one_qr_with_perm(f, x);
    
    f = testclass.make(@(x) [1./(1+1i*x.^2) sinh((1-1i)*x) (exp(x) - x.^3)], ...
        [], [], pref);
    pass(n, 13:14) = test_one_qr(f, x);
    pass(n, 15:16) = test_one_qr_with_perm(f, x);
    
    %%
    % Check that the 'vector' flag works properly.
    N = size(f, 2);
    [Q1, R1, E1] = qr(f);
    [Q2, R2, E2] = qr(f, 'vector');
    err = E1(:, E2) - eye(N);
    pass(n, 17) = all(err(:) == 0);
end

%%
% Restore the RNG state.

rng(rngstate);

end

% Tests the QR decomposition for a CHEBTECH object F using a grid of points X
% in [-1  1] for testing samples.
function result = test_one_qr(f, x)
    N = size(f, 2);
    [Q, R] = qr(f);

    % Check orthogonality.
    ip = innerProduct(Q, Q);
    result(1) = max(max(abs(ip - eye(N)))) < 10*f.epslevel;

    % Check that the factorization is accurate.
    err = Q*R - f;
    result(2) = norm(feval(err, x), 'inf') < 100*f.epslevel;
end

% Same as the previous function but this time uses the QR factorization with
% permutations.
function result = test_one_qr_with_perm(f, x)
    N = size(f, 2);
    [Q, R, E] = qr(f);

    % Check orthogonality.
    ip = innerProduct(Q, Q);
    result(1) = max(max(abs(ip - eye(N)))) < 10*f.epslevel;

    % Check that the factorization is accurate.
    err = Q*R - f*E;
    result(2) = norm(feval(err, x), 'inf') < 100*f.epslevel;
end
