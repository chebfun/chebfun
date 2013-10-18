% Test file for chebtech/bary.m

function pass = test_bary(pref)

% Note. Bary is tested fairly extensively by TEST_FEVAL().

if ( nargin < 1 )
    pref = chebtech.pref();
end

tol = 20*pref.eps;

pass = zeros(2, 2); % Pre-allocate pass matrix.
for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else
        testclass = chebtech2();
    end

    k = 14;
    m = 10;

    x = testclass.chebpts(k);
    f = @(x) sin(x);
    y = linspace(-1, 1, m).';
    fx = f(x);
    fy = f(y);

    % Second kind formula
    pass(n, 1) = norm( testclass.bary(y, fx) - fy ) < tol;
    pass(n, 2) = norm( testclass.bary(y, [fx fx]) - [fy fy] ) < tol;
end

end
