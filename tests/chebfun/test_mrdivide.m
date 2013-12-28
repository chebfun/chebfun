function pass = test_mrdivide(pref)

if ( nargin == 0 )
    pref = chebpref();
end

T = restrict(chebpoly(0:3), [-1 -0.5 0 0.5 1]);
L = restrict(legpoly(0:3), [-1 0 1]);

%% Scalar A:
pass(1) = normest(T/2 - .5*T) < 10*epslevel(T);

%% Numeric A:
B = T;
A = (1:4);
x = A/B;
x0 = feval(x, 0);
x0_true = -2.625;
pass(2) = abs(x0 - x0_true) < 10*epslevel(x);

A = eye(4);
X = A/L;
X0 = feval(X, 0);
pass(3) = norm(X0 - [.5 0 -1.25 0].') < 10*epslevel(X);

%% ~Transposed A:
A = T;
B = (1:4);
x = A/B;
x0 = feval(x, 0);
x0_true = -1/15;
pass(4) = abs(x0 - x0_true) < 10*epslevel(x);

%% Else:
X = T.'/L.';
C = diag([1 1 4/3 8/5]); 
C(3, 1) = -1/3;
C(4,2) = -3/5;
pass(5) = norm(X - C, inf) < 1e2*epslevel(L);

end
