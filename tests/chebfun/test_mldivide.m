function pass = test_mldivide(pref)

if ( nargin == 0 )
    pref = chebfun.pref();
end

T = chebpoly(0:3);
L = legpoly(0:3);

%% Scalar A:
pass(1) = normest(2\T - .5*T) < 10*epslevel(T);

%% Numeric A:
B = T.';
A = (1:4).';
x = A\B;
x0 = feval(x, 0);
x0_true = -1/15;
pass(2) = abs(x0 - x0_true) < 10*epslevel(x);

A = eye(4);
X = A\B;
X0 = feval(X, 0);
pass(3) = norm(X0 - [1 0 -1 0].') < 10*epslevel(X);

%% Transposed A:
A = T.';
B = (1:4).';
X = A\B;
X0 = feval(X, 0);
X0_true = -2.625;
pass(4) = abs(X0 - X0_true) < 10*epslevel(X);

%% Else

X = T\L;
C = diag([1 1 .75 .625]); 
C(1, 3) = .25;
C(2, 4) = .375;
pass(5) = norm(X - C, inf) < 10*epslevel(L);

end
 

