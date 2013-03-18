% Test file for funcheb/bary.

function pass = test_bary(pref)

% [TODO]: This test is rubbish!

if ( nargin < 1 )
    pref = funcheb.pref;
end

tol = 20*pref.funcheb.eps;

for ( n = 1:2 )
    if ( n == 1 )
        testclass = funcheb1();
    else
        testclass = funcheb2();
    end

    k = 14;
    m = 10;

    x = testclass.chebpts(k);
    f = @(x) sin(x);
    y = linspace(-1, 1, m).';
    fx = f(x);
    fy = f(y);

    % second kind formula
    pass(n, 1) = norm( testclass.bary(y, fx, 2) - fy ) < tol;
    pass(n, 2) = norm( testclass.bary(y, [fx fx], 2) - [fy fy] ) < tol;

    % First kind formula
    pass(n, 3) = norm( testclass.bary(y, fx, 1) - fy ) < tol;
    pass(n, 4) = norm( testclass.bary(y, [fx fx], 1) - [fy fy] ) < tol;
end

end
