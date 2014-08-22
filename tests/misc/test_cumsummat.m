function pass = test_cumsummat(pref)

if ( nargin == 0 )
    pref = cheboppref();
end

tol = 1e-10;

c1 = chebcolloc1();
c2 = chebcolloc2();


N = 5;
Da = cumsummat(N);
Db = c2.cumsummat(N);
err = norm(Da - Db);
pass(1) = err < tol;

N = 5;
Da = cumsummat(N, [-2 2]);
Db = c2.cumsummat(N)*2;
err = norm(Da - Db);
pass(2) = err < tol;

N = 5;
Da = cumsummat(N, [-2 2], c2);
Db = c2.cumsummat(N)*2;
err = norm(Da - Db);
pass(3) = err < tol;

N = 5;
Da = cumsummat(N, @chebcolloc1);
Db = c1.cumsummat(N);
err = norm(Da - Db);
pass(4) = err < tol;

N = 5;
Da = cumsummat(N, [-.5 .5], 'chebcolloc1');
Db = c1.cumsummat(N)/2;
err = norm(Da - Db);
pass(5) = err < tol;

p = cheboppref();
cheboppref.setDefaults('discretization', @chebcolloc2);
try
    N = 5;
    Da = cumsummat(N);
    Db = c2.cumsummat(N);
    err = norm(Da - Db);
    pass(6) = err < tol;
catch ME
    chebopref.setDefaults(p);
    rethrow(ME)
end
cheboppref.setDefaults(p);

end
