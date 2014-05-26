function pass = test_splitting(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% Store current globlal state.
p = chebfunpref();
chebfunpref.setDefaults('factory')

try
    % Testing takes place inside this TRY CATCH.
   
    pass(1) = strcmp(splitting(), 'off');
    
    pass(2) = strcmp(splitting('on'), 'off');
    
    pass(3) = strcmp(splitting(), 'on');
    
    F = @(x) cos(pi*x).*sign(x);
    f = chebfun(@(x) F(x), [-1 1]);
    x = linspace(-1, 1, 20);
    err = norm(f(x) - F(x), inf);
    pass(4) = numel(f.domain > 1) && err < 2*epslevel(f);
    
    pass(5) = strcmp(splitting('off'), 'on');
    
    pass(6) = strcmp(splitting(), 'off');
    
catch ME
    chebfunpref.setDefaults(p);
    rethrow(ME)
end

% Return to default settings:
chebfunpref.setDefaults(p)

end