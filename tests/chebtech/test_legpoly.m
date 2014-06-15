% Test file for chebtech/legpoly.m

function pass = test_legpoly(pref)

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
    pass(n, 1) = norm(legpoly(f) - [1 0 0]', inf) < tol;
    
    % Test a vector function (with [P_1, P_2, P_3])
    x = testclass.chebpts(3);
    f = testclass.make([.5*(3*x.^2 - 1), x, (1 + 0*x)], [], pref);
    pass(n, 2) = norm(legpoly(f) - eye(3), inf) < tol;
    
end
