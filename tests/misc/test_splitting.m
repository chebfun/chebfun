function pass = test_splitting(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% Store current globlal state.
p = chebfunpref();
chebfunpref.setDefaults('factory')

warnState = warning('off', 'CHEBFUN:splitting:deprecated');

try
    % Testing takes place inside this TRY CATCH.
   
    pass(1) = strcmp(splitting(), 'off');
    
    pass(2) = strcmp(splitting('on'), 'off');
    
    pass(3) = strcmp(splitting(), 'on');
    
    F = @(x) cos(pi*x).*sign(x);
    f = chebfun(@(x) F(x), [-1 1]);
    x = linspace(-1, 1, 20);
    err = norm(f(x) - F(x), inf);
    pass(4) = numel(f.domain > 1) && err < 1e-5;
    
    
    pass(5) = strcmp(splitting('off'), 'on');
    
    pass(6) = strcmp(splitting(), 'off');
    
catch ME
    % Reset preferences and warning state
    chebfunpref.setDefaults(p);
    warning(warnState);
    
    % Rethrow error:
    rethrow(ME)
end

% Return to default settings:
chebfunpref.setDefaults(p);
warning(warnState);

end
