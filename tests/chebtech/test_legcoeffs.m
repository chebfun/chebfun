% Test file for chebtech/legcoeffs.m

function pass = test_legcoeffs(pref)

if ( nargin < 1 )
    pref = chebtech.techPref();
end

tol = 10*eps;

for n = 1:2
    
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end
    
    % Test a scalar function (with P_3)
    x = testclass.chebpts(3);
    f = testclass.make(.5*(3*x.^2 - 1), [], pref);
    pass(n, 1) = norm(legcoeffs(f) -  [0 0 1]', inf) < tol;
    
    % Test a vector function (with [P_1, P_2, P_3])
    x = testclass.chebpts(3);
    f = testclass.make([(1 + 0*x), x, .5*(3*x.^2 - 1)], [], pref);
    pass(n, 2) = norm(legcoeffs(f) - eye(3), inf) < tol;
    
end
