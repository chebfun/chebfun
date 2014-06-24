function pass = test_diffmat(pref)

if ( nargin == 0 )
    pref = cheboppref();
end

tol = 1e-10;

c1 = colloc1();
c2 = colloc2();

N = 5;
Da = diffmat(N);
Db = colloc2.diffmat(N);
err = norm(Da - Db);
pass(1) = err < tol;

N = 5;
Da = diffmat(N, [-2 2]);
Db = c2.diffmat(N)/2;
err = norm(Da - Db);
pass(2) = err < tol;

N = 5;
Da = diffmat(N, 2, [-2 2]);
Db = c2.diffmat(N,2)/4;
err = norm(Da - Db);
pass(3) = err < tol;

N = 5;
Da = diffmat(N, c1);
Db = c1.diffmat(N,1);
err = norm(Da - Db);
pass(4) = err < tol;

N = 5;
Da = diffmat(N, 3, [-.5 .5], 'colloc1');
Db = c1.diffmat(N,3)*8;
err = norm(Da - Db);
pass(5) = err < tol;

p = cheboppref();
cheboppref.setDefaults('discretization', @colloc1);
try
    N = 5;
    Da = diffmat(N, 3);
    Db = c2.diffmat(N,3);
    err = norm(Da - Db);
    pass(6) = err < tol;
catch ME
    chebopref.setDefaults(p);
    rethrow(ME)
end
cheboppref.setDefaults(p);

end