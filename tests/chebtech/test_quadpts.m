function pass = test_quadpts(pref)

if ( nargin == 1 )
    pref = chebtech.techPref();
end

for m = 1:2
    if ( m == 1 )
        testTech = chebtech1();
    else
        testTech = chebtech2();
    end
    
    n = 10;
    w = testTech.quadwts(n);
    
    % Test against some low order polynomials:
    x = testTech.chebpts(n);
    pass(m, 1) = abs(sum(w) - 2) < 2*eps;
    pass(m, 2) = abs(w*x) < eps;
    pass(m, 3) = abs(w*x.^2 - 2/3) < 2*eps;
    pass(m, 4) = abs(w*x.^3) < eps;
    pass(m, 5) = abs(w*x.^4 - 2/5) < 2*eps;
    
    % Test against some basic properties:
    pass(m, 6) = isempty(testTech.quadwts(0));
    pass(m, 7) = testTech.quadwts(1) == 2;
    pass(m, 8) = all(testTech.quadwts(2) == 1);
    pass(m, 9) = mod(n,2) || norm(w(1:n/2) - w(n:-1:n/2+1), inf) < eps(n);
    
end

    
