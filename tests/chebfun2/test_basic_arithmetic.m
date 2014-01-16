function pass = test_basic_arithmetic
% This tests the basic arithmetic operations on chebfun2 objects.

tol = 10000 * eps; j = 1;

D = [-1 1 -1 1; -2 2 -2 2; -1 pi 0 2*pi];

for r = 1 : size(D,1)
    f = chebfun2(@(x,y) cos(x.*y), D(r,:));
    g = chebfun2(@(x,y) x + y + x.*y, D(r,:));
    
    uplusF = f;
    uminusF = chebfun2(@(x,y) -cos(x.*y), D(r,:));
    FplusG = chebfun2(@(x,y) cos(x.*y) + (x + y + x.*y), D(r,:) );
    FminusG = chebfun2(@(x,y) cos(x.*y) - (x + y + x.*y), D(r,:) );
    FtimesG = chebfun2(@(x,y) cos(x.*y).*(x + y + x.*y), D(r,:) );
    
    tolr = norm(D(r,:),inf)*tol;
    
    pass(j) = ( norm( f - uplusF ) < tolr ); j = j + 1;
    pass(j) = ( norm( (-f) - uminusF ) < tolr ); j = j + 1;
    pass(j) = ( norm( f+g - FplusG ) < tolr ); j = j + 1;
    pass(j) = ( norm( (f-g) - FminusG ) < 10*tolr ); j = j + 1;
    pass(j) = ( norm( f.*g - FtimesG ) < 10*tolr ); j = j + 1;
end


end