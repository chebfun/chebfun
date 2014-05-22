function pass = test_blowup(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% Store current globlal state.
p = chebfunpref();
chebfunpref.setDefaults('factory')

try
    % Testing takes place inside this TRY CATCH.
   
    pass(1) = blowup() == 0;
    
    pass(2) = blowup('on') == 0;
    
    pass(3) = blowup() == 1;
    
    F = @(x) 1./(x + 1);
    f = chebfun(@(x) F(x) , [-1 1]);
    x = .9*linspace(-1, 1, 20);
    err = norm(f(x) - F(x), inf);
    tol = 1e2*epslevel(f).*vscale(f);
    pass(4) = err < tol;
    
    p2 = chebfunpref();
    pass(5) = all(strcmp(p2.singPrefs.singType, 'pole'));
    
    pass(6) = blowup(2) == 1;
    
    pass(7) = blowup() == 2;
    
    F = @(x) 1./sqrt(x + 1);
    f = chebfun(@(x) F(x) , [-1 1]);
    x = .9*linspace(-1, 1, 20);
    err = norm(f(x) - F(x), inf);
    tol = 1e2*epslevel(f).*vscale(f);
    pass(8) = err < tol;
    
    pass(9) = blowup('off') == 2;
    
    pass(10) = blowup() == 0;
    
catch ME
    chebfunpref.setDefaults(p);
    rethrow(ME)
end

% Return to default settings:
chebfunpref.setDefaults(p)

end