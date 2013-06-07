% Test file for bndfun/qr.m

function pass = test_qr(pref)

% Get preferences.
if ( nargin < 1 )
    pref = bndfun.pref;
end

% Create seed for random number generator
seedRNG(6178);

pass = zeros(2, 17); % Pre-allocate pass matrix
for n = 1:1 %[TODO]: unbndfun
    if ( n == 1 )
        testclass = bndfun();
        dom = [-2 7];
        
        % Generate a few random points to use as test values.
        x = diff(dom) * rand(100, 1) + dom(1);
    else
        testclass = unbndfun();
    end

    %%
    % Do a few spot-checks.
    
    f = testclass.make(@(x) sin(x), dom, [], [], pref);
    pass(n, 1:2) = test_one_qr(f, x);
    pass(n, 3:4) = test_one_qr_with_perm(f, x);
    
    f = testclass.make(@(x) [cos(x) exp(x)], dom, [], [], pref);
    pass(n, 5:6) = test_one_qr(f, x);
    pass(n, 7:8) = test_one_qr_with_perm(f, x);
    
    f = testclass.make(@(x) [ones(size(x)) x x.^2 x.^3 x.^4 x.^5 x.^6 x.^7], ...
        dom, [], [], pref);
    pass(n, 9:10) = test_one_qr(f, x);
    pass(n, 11:12) = test_one_qr_with_perm(f, x);
    
    f = testclass.make(@(x) [1./(1+1i*x.^2) sinh((1-1i)*x) (exp(x) - x.^3)], ...
        dom, [], [], pref);
    pass(n, 13:14) = test_one_qr(f, x);
    pass(n, 15:16) = test_one_qr_with_perm(f, x);
    
    %%
    % Check that the 'vector' flag works properly.
    N = size(f, 2);
    [Q1, R1, E1] = qr(f);
    [Q2, R2, E2] = qr(f, 'vector');
    err = E1(:, E2) - eye(N);
    pass(n, 17) = all(err(:) == 0);
    
    %%
    % Check a rank-deficient problem:
    f = testclass.make(@(x) [x x x], dom, [], [], pref);
    [Q, R] = qr(f);
    pass(n, 18) = all(size(Q) == 3) && all(size(R) == 3);
    I = eye(3); I(end) = 0;
    pass(n, 19) = norm(innerProduct(Q, Q) - I, inf) < 10*f.epslevel;
    
end

end

% Tests the QR decomposition for a CHEBTECH object F using a grid of points X
% in [-1  1] for testing samples.
function result = test_one_qr(f, x)
    N = size(f, 2);
    [Q, R] = qr(f);

    % Check orthogonality.
    ip = innerProduct(Q, Q);
    result(1) = max(max(abs(ip - eye(N)))) < max(f.onefun.vscale)*f.onefun.epslevel;

    % Check that the factorization is accurate.
    err = Q*R - f;
    result(2) = norm(feval(err, x), inf) < max(f.onefun.vscale)*f.onefun.epslevel;
end

% Same as the previous function but this time uses the QR factorization with
% permutations.
function result = test_one_qr_with_perm(f, x)
    N = size(f, 2);
    [Q, R, E] = qr(f);

    % Check orthogonality.
    ip = innerProduct(Q, Q);
    result(1) = max(max(abs(ip - eye(N)))) < max(f.onefun.vscale)*f.onefun.epslevel;

    % Check that the factorization is accurate.
    err = Q*R - f*E;
    result(2) = norm(feval(err, x), inf) < max(f.onefun.vscale)*f.onefun.epslevel;
end
