function pass = test_vectorCheck(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

warnState = warning();
warning('off', 'CHEBFUN:CHEBFUN:vectorcheck:vectorize');
warning('off', 'CHEBFUN:CHEBFUN:vectorcheck:transpose');
    
try 

    f = chebfun(@(x) 1, pref);
    g = chebfun(@(x) 1+0*x, pref);
    pass(1) = norm(f - g) == 0;

    f = chebfun(@(x) sin(x).', pref);
    g = chebfun(@sin, pref);
    pass(2) = norm(f - g) == 0;

    f = chebfun(@(x) sin(x).');
    pass(2) = norm(f - chebfun(@sin)) == 0;

    f = chebfun(@(x) quadgk(@cos, x, 2));
    g = chebfun(@(x) quadgk(@cos, x, 2), 'vectorize');
    pass(3) = norm(f - g) == 0;
    
    f = chebfun(@(x) [1 2 3], pref);
    g = chebfun([1 2 3], pref);
    pass(4) = norm(f - g) == 0;
    
    f = chebfun(@(x) [x x].', [-1 1], pref);
    g = chebfun(@(x) [x x], [-1 1], pref);
    pass(5) = norm(f - g) == 0;
    
    f = chebfun(@(x) [x x].', [-1 0 1], pref);
    g = chebfun(@(x) [x x], [-1 0 1], pref);
    pass(6) = norm(f - g) == 0;
    
    f = chebfun(@(x) [x x x].', [-1 1], pref);
    g = chebfun(@(x) [x x x], [-1 1], pref);
    pass(7) = norm(f - g) == 0;
    
    f = chebfun(@(x) [exp(-x.^2) exp(-x.^2)].', [0, inf], pref);
    g = chebfun(@(x) [exp(-x.^2) exp(-x.^2)], [0, inf], pref);
    pass(8) = normest(f - g) < 1e-10;
    
    % Tests from #937
    f = chebfun(@(x)[sin(x) 1], pref);
    g = chebfun(@(x) [sin(x) 1+0*x], pref);
    pass(9) = norm(f - g) == 0;
    
    f = chebfun(@(x)[1 1], pref);
    g = chebfun([1 1], pref);
    pass(10) = norm(f - g) == 0;
    
    % See #1671.
    f = chebfun(@(x) x/cos(x), pref);
    g = chebfun(@(x) x./cos(x), pref);
    pass(11) = norm(f - g) == 0;
    
    warning(warnState);
    
catch ME
    
    warning(warnState);
    rethrow(ME);
    
end
    
end
