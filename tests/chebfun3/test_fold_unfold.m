function pass = test_fold_unfold()
% Test chebfun3/fold and chebfun3/unfold

m = 5; n = 4; p = 6;
T = rand(m, n, p);

% Mode 1:
M1 = chebfun3.unfold(T, 1); % Mode-1 unfolding of T
pass(1) = all(size(M1) == [m n*p]);

% Check the reverse operation:
T_new = chebfun3.fold(M1, [m, n, p], 1, [2 3]);
err = T_new - T;
pass(2) = norm(err(:)) == 0;

% Mode 2:
M2 = chebfun3.unfold(T, 2); % Mode-2 unfolding of T
pass(3) = all(size(M2) == [n m*p]);

% Check the reverse operation:
T_new = chebfun3.fold(M2, [m, n, p], 2, [1 3]);
err = T_new - T;
pass(4) = norm(err(:)) == 0;

% Mode 3:
M3 = chebfun3.unfold(T, 3); % Mode-3 unfolding of T
pass(5) = all(size(M3) == [p m*n]);

% Check the reverse operation:
T_new = chebfun3.fold(M3, [m, n, p], 3, [1 2]);
err = T_new - T;
pass(6) = norm(err(:)) == 0;

end